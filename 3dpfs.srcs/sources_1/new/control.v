`timescale 1ns / 1ps

module control #(
	parameter LEN_BITS = 8,
	parameter RECV_BUF_BITS = 10,
	parameter SPI_DIVIDER = 100,
	parameter SPI_BITS_PER_CHIP = 40,
	parameter SPI_CHIPS = 3,
	parameter SPIBITS = SPI_BITS_PER_CHIP * SPI_CHIPS,
	parameter REGBITS = SPIBITS,
	parameter NCS = 2,
	parameter NGPOUT = 1,
	parameter NGPIN = 0
) (
	input clk,

	/* len fifo input */
	input len_fifo_empty,
	input [LEN_BITS-1:0] len_fifo_data,
	output reg len_fifo_rd_en,

	/* ring buffer input */
	input [7:0] recv_data,	/* data at rptr */
	output reg [RECV_BUF_BITS-1:0] recv_rptr = 0,

	/* stepper board control */
	output reg sck = 0,
	output reg [NCS-1:0] cs = { (NCS){ 1'b1 } },
	input sdo,		/* output from slave */
	output reg sdi = 0,	/* input to slave */
	output reg [NGPOUT-1:0] gpout = 0,
	input [NGPIN-1:0] gpin,

	/* debug output */
	output [31:0] debug
);

parameter SPI_DIVIDER_BITS = $clog2(SPI_DIVIDER);
parameter SPICNT_WIDTH = $clog2(SPIBITS);
parameter NGPOUT_BITS = $clog2(NGPOUT);
parameter NGPIN_BITS = $clog2(NGPIN);

localparam C_CMD_SPI = 8'h80;	/* 0x80-0x8f, lower 4 bit for cs */
localparam C_CMD_GPOUT_HI = 8'h70;
localparam C_CMD_GPOUT_LO = 8'h71;
localparam C_CMD_GPIN = 8'h78;

reg [7:0] c_len = 0;
reg [7:0] c_cmd = 0;
reg [REGBITS-1:0] c_creg = 0;
reg c_in_cmd = 0;
reg c_cmd_end = 0;

reg c_spi_start = 0;
reg s_spi_running = 0;
reg s_spi_done = 0;

always @(posedge clk) begin: main_block
	integer i;

	len_fifo_rd_en <= 0;
	if (c_len == 0 && !len_fifo_empty && !c_in_cmd) begin
		/*
		 * stage: read command byte
		 */
		c_cmd <= recv_data;
		c_in_cmd <= 1;
		c_len <= len_fifo_data - 1;
		c_creg <= 0;
		recv_rptr <= recv_rptr + 1;
	end else if (c_len != 0) begin
		/*
		 * stage: compose value
		 */
		c_creg <= { c_creg[REGBITS-9:0], recv_data };
		c_len <= c_len - 1;
		recv_rptr <= recv_rptr + 1;
	end else if (c_in_cmd && !c_cmd_end) begin
		/*
		 * last stage: interpret command
		 */
		if (c_spi_start || s_spi_running) begin
			/* do nothing while spi is running */
			c_spi_start <= 0;
		end else if (s_spi_done) begin
			/* spi finished */
			c_cmd_end <= 1;
			len_fifo_rd_en <= 1;
			/* TODO: send result */
		end else if (c_cmd[7:4] == C_CMD_SPI[7:4]) begin
			c_spi_start <= 1;
		end else if (c_cmd == C_CMD_GPOUT_HI) begin
			gpout[c_creg[NGPOUT-1:0]] <= 1;
			c_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (c_cmd == C_CMD_GPOUT_LO) begin
			gpout[c_creg[NGPOUT-1:0]] <= 0;
			c_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else begin
			/*
			 * unknwon command, ignore
			 */
		end
	end else if (c_cmd_end) begin
		/*
		 * end of command stage. delayed by one clock so the recv
		 * len fifo can advance
		 */
		c_cmd_end <= 0;
		c_in_cmd <= 0;
	end
end

/*
 * SPI interface
 */
reg [SPI_DIVIDER_BITS-1:0] s_spi_divider = 0;
reg [SPIBITS-1:0] s_spi_data_out = 0;
reg [SPIBITS-1:0] s_spi_data_in = 0;
reg [SPICNT_WIDTH-1:0] s_spi_cnt = 0;
reg [1:0] s_spi_phase = 0;
reg sdo_sync1 = 0;
reg sdo_sync = 0;
always @(posedge clk) begin
	sdo_sync1 <= sdo;
	sdo_sync <= sdo_sync1;
end

always @(posedge clk) begin
	if (c_spi_start && !s_spi_running) begin
		sck <= 1;
		s_spi_cnt <= SPIBITS;
		s_spi_running <= 1;
		s_spi_phase <= 0;
		s_spi_divider <= SPI_DIVIDER;
		s_spi_data_out <= c_creg[SPIBITS-1:0];
		s_spi_data_in <= 0;
	end else if (s_spi_divider != 0) begin
		s_spi_divider <= s_spi_divider - 1;
	end else if (s_spi_running && s_spi_phase == 0) begin
		/* phase 0: assert CS */
		cs[1 << c_cmd[3:0]] <= 0;
		s_spi_divider <= SPI_DIVIDER;
		s_spi_phase <= 1;
	end else if (s_spi_running && s_spi_phase == 1) begin
		/* phase 1: falling clock, output bit */
		sck <= 0;
		sdi <= s_spi_data_out[SPIBITS-1];
		s_spi_data_out <= { s_spi_data_out[SPIBITS-2:0], 1'h0 };
		s_spi_divider <= SPI_DIVIDER;
		s_spi_phase <= 2;
	end else if (s_spi_running && s_spi_phase == 2) begin
		/* phase 2: rising clock, sample input */
		sck <= 1;
		s_spi_data_in <= { s_spi_data_in[SPIBITS-2:0], sdo_sync };
		s_spi_divider <= SPI_DIVIDER;
		s_spi_cnt <= s_spi_cnt - 1;
		if (s_spi_cnt == 1)
			s_spi_phase <= 3;
		else
			s_spi_phase <= 1;
	end else if (s_spi_running && s_spi_phase == 3) begin
		cs[1 << c_cmd[3:0]] <= 1;
		s_spi_done <= 1;
		s_spi_running <= 0;
	end else if (s_spi_done) begin
		s_spi_done <= 0;
	end 
end

endmodule
