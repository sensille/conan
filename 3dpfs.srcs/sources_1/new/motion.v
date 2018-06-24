`timescale 1ns / 1ps

module motion #(
	parameter LEN_BITS = 8,
	parameter RECV_BUF_BITS = 10,
	parameter REGBITS = 64,
	parameter NCNTRL = 4
) (
	input clk,

	/* len fifo input */
	input len_fifo_empty,
	input [LEN_BITS-1:0] len_fifo_data,
	output reg len_fifo_rd_en,

	/* ring buffer input */
	input [7:0] recv_data,	/* data at rptr */
	output reg [RECV_BUF_BITS-1:0] recv_rptr = 0,

	/* step/dir output */
	output reg [NCNTRL-1:0] step = 0,
	output [NCNTRL-1:0] dir,

	/* debug output */
	output [31:0] debug
);

/*
 * Preload registers for motion stage
 */
reg [REGBITS-1:0] m_preload_reg [NCNTRL-1:0];

/*
 * motion state machine
 */
localparam M_STATE_STOP = 0;
localparam M_STATE_RUN = 1;
localparam M_STATE_FREERUN = 2;
localparam M_STATE_ERROR = 3;
reg [1:0] m_state = M_STATE_STOP;
reg [REGBITS-1:0] m_cnt = 0;
reg [REGBITS-1:0] m_jerk [NCNTRL-1:0];
reg [REGBITS-1:0] m_accel [NCNTRL-1:0];
reg [REGBITS-1:0] m_velocity [NCNTRL-1:0];
reg [REGBITS-1:0] m_pos [NCNTRL-1:0];

/* initialize array to 0, mainly for simulation */
initial begin: init_mem
	integer i;
	for (i = 0; i < NCNTRL - 1; i = i + 1) begin
		m_jerk[i] = 0;
		m_accel[i] = 0;
		m_velocity[i] = 0;
		m_pos[i] = 0;
	end
end
/*
 * packet parsing state machine
 */
localparam P_STATE_IDLE = 0;
localparam P_STATE_BUILD_VAL = 1;

/*
 * command definitions
 */
localparam M_CMD_LOADALLREG = 8'h70;
localparam M_CMD_LOADCNT    = 8'h71;
localparam M_CMD_LOADREG    = 8'h80;	/* 0x80-0x8f */
/*
 * Queue handling, read packet from queue, interpret command,
 * fill preload registers
 */
reg [LEN_BITS-1:0] mp_len = 0;	/* reading len bytes */
reg [REGBITS-1:0] mp_creg = 0;	/* composition register, build value */
reg [7:0] mp_cmd = 0;		/* current cmd */
reg mp_in_cmd = 0;		/* currently parsing cmd */
reg mp_cmd_end = 0;		/* end of command delay clk */
reg mp_extend_sign = 0;		/* set when first byte is received */
reg i_want_an_empty_block;	/* keep compiler happy */
always @(posedge clk) begin
	len_fifo_rd_en <= 0;
	if (m_state == M_STATE_ERROR) begin
		/*
		 * just stay in error state
		 */
		i_want_an_empty_block <= 0;
	end else if (mp_len == 0 && !len_fifo_empty && !mp_in_cmd) begin
		/*
		 * stage: read command byte
		 */
		mp_cmd <= recv_data;
		mp_in_cmd <= 1;
		mp_len <= len_fifo_data - 1;
		mp_creg <= 0;
		mp_extend_sign <= 1;
		recv_rptr <= recv_rptr + 1;
	end else if (mp_len != 0) begin
		/*
		 * stage: compose value
		 */
		if (mp_extend_sign) 	/* first byte, extend sign */
			mp_creg <= { {(REGBITS-8){recv_data[7]}}, recv_data };
		else			/* all other bytes, shift data in */
			mp_creg <= { mp_creg[REGBITS-9:0], recv_data };
		mp_len <= mp_len - 1;
		mp_extend_sign <= 0;
		recv_rptr <= recv_rptr + 1;
	end else if (mp_in_cmd && !mp_cmd_end) begin
		/*
		 * last stage: interpret command
		 */
		if (mp_cmd[7:4] == M_CMD_LOADREG[7:4]) begin
			m_preload_reg[mp_cmd[3:0]] = mp_creg;
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_LOADALLREG) begin: loadall
			integer i;
			for (i = 0; i < NCNTRL; i = i + 1) begin
				m_preload_reg[i] <= mp_creg;
			end
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_LOADCNT) begin
			if (m_state == M_STATE_RUN && m_cnt != 0) begin
				/*
				 * current motion still running, wait for it
				 * to finish before loading cnt
				 */
				i_want_an_empty_block <= 0;
			end else if (m_state == M_STATE_FREERUN) begin
				/*
				 * current motion runs until an event occurs
				 */
				/* XXX TODO */
				i_want_an_empty_block <= 0;
			end else if (m_state != M_STATE_ERROR) begin: loadcnt
				integer i;
				/*
				 * load staged registers + cnt
				 */
				m_cnt <= mp_creg;
				for (i = 0; i < NCNTRL; i = i + 1) begin
					m_jerk[i] <= m_preload_reg[i];
				end;
				mp_cmd_end <= 1;
				len_fifo_rd_en <= 1;
				m_state <= M_STATE_RUN;
			end
		end else begin
			/*
			 * unknwon command
			 */
			m_state <= M_STATE_ERROR;
		end
			
	end else if (mp_cmd_end) begin
		/*
		 * end of command stage. delayed by one clock so the recv
		 * len fifo can advance
		 */
		mp_cmd_end <= 0;
		mp_in_cmd <= 0;
	end
end

/*
 * motion controller
 */
reg [NCNTRL-1:0] do_step;
genvar gi;
generate
	for (gi = 0; gi < NCNTRL; gi = gi + 1) begin
		assign dir[gi] = m_velocity[gi][REGBITS-1];
	end
endgenerate 
always @(posedge clk) begin: motion_block
	integer i;

	for (i = 0; i < NCNTRL; i = i + 1) begin
		if (do_step[i]) begin
			step[i] <= ~step[i];
		end
		do_step[i] <= 0;
	end
	if (m_cnt != 0 || m_state == M_STATE_FREERUN) begin
		for (i = 0; i < NCNTRL; i = i + 1) begin
			m_accel[i] <= m_accel[i] + m_jerk[i];
			m_velocity[i] <= m_velocity[i] + m_accel[i];
			/* do_step is the overflow of the addition */
			{ do_step[i], m_pos[i] } <=
				{ m_pos[i][REGBITS-1], m_pos[i] } +
			        { m_velocity[i][REGBITS-1], m_velocity[i] };
		end
		m_cnt <= m_cnt - 1;
	end
end

assign debug = { m_state, mp_cmd, mp_in_cmd };

endmodule
