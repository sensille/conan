`timescale 1ns / 1ps

module main(
	input clk_50mhz,

	output led,

	input rx,
	output tx,

	// stepper interface
	output en,
	output sck,
	output cs123,
	output cs456,
	output sdi,
	input sdo,
	output [1:6] dir,
	output [1:6] step,

	// debug
	output debug1,
	output debug2,
	output debug3,
	output debug4,

	//
	output leds_out,
	output leds_clk,
	output leds_cs
);

localparam NCNTRL = 6;
localparam JERK_WIDTH = 36;
localparam CNT_WIDTH = 36;
localparam RECV_BUF_BITS = 10;
localparam RECV_BUF_SIZE = (1 << 10);

assign debug1 = rx;
assign debug2 = tx;

wire clk;
clk u_clk(
	.clk_50mhz(clk_50mhz),
	.clk(clk)		// 20 MHz
);

reg [21:0] led_cnt = 0;
always @(posedge clk)
	led_cnt = led_cnt + 1;
assign led = led_cnt[21];

reg rx_s1 = 0;
reg rx_sync = 0;
always @(posedge clk) begin
	rx_s1 = rx;
	rx_sync = rx_s1;
end

wire [7:0] tx_data;
wire [7:0] rx_data;
wire tx_en;
wire rx_ready;

parameter CLOCK_DIVIDE = 521; // clock rate (20Mhz) / (baud rate (9600) * 4)
uart #(
        .CLOCK_DIVIDE(CLOCK_DIVIDE)
) uart_u (
	.clk(clk),
	.rst(0),
	.rx(rx_sync),
	.tx(tx),
	.transmit(tx_en),
	.tx_byte(tx_data),
	.received(rx_ready),
	.rx_byte(rx_data),
	.is_receiving(),
	.is_transmitting(),
	.recv_error()
);

fifo #(
	.DATA_WIDTH(8)
) recv_data_fifo (
	.clk(clk),
	.clr(0),

	// write side
	.din(rx_data),
	.wr_en(fifo_wr_en),
	.full(),

	// read side
	.dout(fifo_rd_data),
	.rd_en(fifo_rd_en_combined),
	.empty(fifo_rd_empty),

	// status
	.elemcnt()
);

reg [7:0] recv_ring [RECV_BUF_SIZE-1:0];
reg [RECV_BUF_BITS-1:0] recv_wptr = 0;
reg [RECV_BUF_BITS-1:0] recv_wptr_fallback = 0;
reg [RECV_BUF_BITS-1:0] recv_rptr = 0;
reg recv_state = 0;

localparam RECV_IDLE = 0;
localparam RECV_IN_FRAME = 1;
localparam RECV_IN_ESCAPE = 2;
localparam RX_ESC_CHAR = 8'h7d;

wire recv_data = (recv_state == RECV_IN_ESCAPE) ? rx_data ^ 8'h20 : rx_data;
wire next_recv_wptr = recv_wptr + 1;	/* wraps */

always @(posedge clk) begin
	if (rx_ready) begin
		if (recv_state == RECV_IDLE) begin
			recv_wptr_fallback = recv_wptr;
		end
		if (rx_data == RECV_ESC_CHAR) begin
			recv_state = RECV_IN_ESCAPE;
		end else if (rx_data == RECV_EOF_CHAR) begin
			if (recv_state == RECV_IN_ESCAPE ||
			    calc_crc != recv_crc) begin
				rx_error = 1;
				recv_wptr = recv_wptr_fallback;
			end else begin
				/* push len */
			end
		end else begin
			recv_ring[recv_wptr] = recv_data;
			recv_wptr = next_recv_wptr;
		end
	end
end

assign tx_data = rx_data ^ 8'h23;
assign tx_en = rx_ready;

wire di = 0;
wire crc0 = di ^ crc16[15];
reg [15:0] crc16;
always @(posedge clk) begin
	crc16[0] = crc0;
	crc16[1] = crc16[0];
	crc16[2] = crc16[1];
	crc16[3] = crc16[2];
	crc16[4] = crc16[3];
	crc16[5] = crc16[4] ^ crc0;
	crc16[6] = crc16[5];
	crc16[7] = crc16[6];
	crc16[8] = crc16[7];
	crc16[9] = crc16[8];
	crc16[10] = crc16[9];
	crc16[11] = crc16[10];
	crc16[12] = crc16[11] ^ crc0;
	crc16[13] = crc16[12];
	crc16[14] = crc16[13];
	crc16[15] = crc16[14];
end

wire abort = 0;
wire [NCNTRL:1] fifo_rd_en;
wire fio_rd_en_combined = |fifo_rd_en;
wire [JERK_WIDTH + CNT_WIDTH + NCNTRL - 1 : 0] fifo_wr_data;
wire [JERK_WIDTH + CNT_WIDTH + NCNTRL - 1 : 0] fifo_rd_data;
wire fifo_rd_empty;

wire [6:0] elemcnt;
fifo #(
	.DATA_WIDTH(JERK_WIDTH + CNT_WIDTH + NCNTRL)
) u_fifo (
	.clk(clk),
	.clr(abort),

	// write side
	.din(fifo_wr_data),
	.wr_en(fifo_wr_en),
	.full(),

	// read side
	.dout(fifo_rd_data),
	.rd_en(fifo_rd_en_combined),
	.empty(fifo_rd_empty),

	// status
	.elemcnt()
);


// TODO: check channel

/*
genvar k;
generate for (k = 1; k < NCNTRL; k = k + 1) begin
	cntrl cntrl_u(
		.clk(clk),
		.rd_data(fifo_rd_data),
		.rd_empty(fifo_rd_empty),
		.rd_en(fifo_rd_en[k]),
		.step(step[k]),
		.dir(dir[k]),
		.start(start),
		.abort(abort),
		.underflow(),
		.running(),
		.idle(),
		.current_velocity()
	);
end
*/

reg [255:0] leds = 256'h123456789abcdef023456789abcdef013456789abcdef012456789abcdef0123;
led7219 led7219_u(
	.clk(clk),
	.data(leds),
	.leds_out(leds_out),
	.leds_clk(leds_clk),
	.leds_cs(leds_cs)
);

endmodule
