`timescale 1ns / 1ps
`default_nettype none

module motion #(
	parameter LEN_BITS = 8,
	parameter RECV_BUF_BITS = 10,
	parameter REGBITS = 208,
	parameter CNTBITS = 48,
	parameter NCNTRL = 4,
	parameter NSTEPDIR = 6,
	parameter NENDSTOP = 2
) (
	input wire clk,

	/* len fifo input */
	input wire len_fifo_empty,
	input wire [LEN_BITS-1:0] len_fifo_data,
	output reg len_fifo_rd_en,

	/* ring buffer input */
	input wire [7:0] recv_data,	/* data at rptr */
	output reg [RECV_BUF_BITS-1:0] recv_rptr = 0,

	/* send len fifo */
	output reg send_fifo_wr_en = 0,
	output reg [LEN_BITS-1:0] send_fifo_data = 0,
	input wire send_fifo_full,

	/* send ring */
	output reg [7:0] send_ring_data = 0,
	output reg send_ring_wr_en = 0,
	input wire send_ring_full,

	/* step/dir output */
	output reg [NSTEPDIR-1:0] step,
	output reg [NSTEPDIR-1:0] dir,

	/* run control */
	input wire running,
	input wire clear_error,
	output reg request_stop = 0,

	/* endstop sensors */
	input wire [NENDSTOP-1:0] endstop,

	/* debug output */
	output wire [31:0] debug
);
localparam NCNTRL_BITS = $clog2(NCNTRL);
localparam NCNTRL_BITS_R = $clog2(NCNTRL + 1); /* +1 for the routing to NULL */
localparam NSTEPDIR_BITS = $clog2(NSTEPDIR);
/* step/dir from controllers */
reg [NCNTRL:0] c_step = 0;
wire [NCNTRL:0] c_dir;
assign c_dir[0] = 0;

/* routing information */
reg [NCNTRL_BITS_R-1:0] stepdir_routing[NSTEPDIR-1:0];
/* initialize array to 0, mainly for simulation */
initial begin: init_routing
	integer i;
	for (i = 0; i < NSTEPDIR; i = i + 1) begin
		stepdir_routing[i] <= 0;
	end
end

/*
 * Preload registers for motion stage
 */
reg [REGBITS-1:0] m_preload_jerk [NCNTRL-1:0];
reg [REGBITS-1:0] m_preload_snap [NCNTRL-1:0];
reg [REGBITS-1:0] m_preload_crackle [NCNTRL-1:0];

/*
 * motion state machine
 */
reg m_state_error = 0;
reg [CNTBITS-1:0] m_cnt = 0;
reg [REGBITS-1:0] m_crackle [NCNTRL-1:0];
reg [REGBITS-1:0] m_snap [NCNTRL-1:0];
reg [REGBITS-1:0] m_jerk [NCNTRL-1:0];
reg [REGBITS-1:0] m_accel [NCNTRL-1:0];
reg [REGBITS-1:0] m_velocity [NCNTRL-1:0];
reg [REGBITS-1:0] m_pos [NCNTRL-1:0];

/* initialize array to 0, mainly for simulation */
initial begin: init_mem
	integer i;
	for (i = 0; i < NCNTRL; i = i + 1) begin
		m_crackle[i] <= 0;
		m_snap[i] <= 0;
		m_jerk[i] <= 0;
		m_accel[i] <= 0;
		m_velocity[i] <= 0;
		m_pos[i] <= 0;
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
localparam M_CMD_RESET      = 8'h40;
localparam M_CMD_STOP       = 8'h41;
localparam M_CMD_SET_ROUTING= 8'h60;
localparam M_CMD_NOTIFY     = 8'h61;
localparam M_CMD_WAIT_IDLE  = 8'h68;
localparam M_CMD_SET_EVENT  = 8'h69;
localparam M_CMD_ENDSTOP_MASK=8'h6a;
localparam M_CMD_ENDSTOP_POL= 8'h6b;
localparam M_CMD_LOADCNT    = 8'h71;
localparam M_CMD_LOADJERK   = 8'h80;	/* 0x80-0x8f */
localparam M_CMD_LOADSNAP   = 8'h90;	/* 0x90-0x9f */
localparam M_CMD_LOADCRACKLE= 8'ha0;	/* 0xa0-0xaf */
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
reg [NCNTRL-1:0] do_step;
reg send_notify = 0;
reg [31:0] send_notify_data = 0;
reg motion_reset = 0;
reg [NENDSTOP-1:0] endstop_mask = 0;
reg [NENDSTOP-1:0] endstop_pol = 0;
wire [NENDSTOP-1:0] endstop_int = endstop ^ endstop_pol;
reg [NENDSTOP-1:0] m_wait_mask = 0;

always @(posedge clk) begin: main_block
	integer i;

	request_stop <= 0;
	len_fifo_rd_en <= 0;
	if (m_cnt != 0)
		m_cnt <= m_cnt - 1;
	else if (running && !request_stop)
		/* underrun, overridden by load if in same cycle */
		m_state_error <= 1;

	if (send_notify)
		send_notify <= 0;

	if (clear_error && m_state_error) begin
		m_state_error <= 0;
		motion_reset <= 1;
		mp_cmd_end <= 0;
		mp_in_cmd <= 0;
		recv_rptr <= 0;
		m_wait_mask <= 0;
	end

	if (m_state_error) begin
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
		if (mp_cmd[7:4] == M_CMD_LOADJERK[7:4]) begin
			m_preload_jerk[mp_cmd[NCNTRL_BITS-1:0]] <= mp_creg;
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd[7:4] == M_CMD_LOADSNAP[7:4]) begin
			m_preload_snap[mp_cmd[NCNTRL_BITS-1:0]] <= mp_creg;
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd[7:4] == M_CMD_LOADCRACKLE[7:4]) begin
			m_preload_crackle[mp_cmd[NCNTRL_BITS-1:0]] <= mp_creg;
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_LOADCNT) begin
			if (!running) begin
				/* wait for running before loading cnt */
				i_want_an_empty_block <= 0;
			end else if (m_cnt != 0 &&
			    (m_wait_mask & endstop_int) == 0) begin
				/*
				 * current motion still running, wait for it
				 * to finish or an event to occur before
				 * loading cnt
				 */
				i_want_an_empty_block <= 0;
			end else begin
				/*
				 * load staged registers + cnt
				 */
				m_cnt <= mp_creg;
				for (i = 0; i < NCNTRL; i = i + 1) begin
					m_jerk[i] <= m_preload_jerk[i];
					m_snap[i] <= m_preload_snap[i];
					m_crackle[i] <= m_preload_crackle[i];
				end;
				mp_cmd_end <= 1;
				len_fifo_rd_en <= 1;
				m_state_error <= 0;
				m_wait_mask <= 0;
			end
		end else if (mp_cmd == M_CMD_SET_ROUTING) begin
			stepdir_routing[mp_creg[NSTEPDIR_BITS-1+8:8]] <=
				mp_creg[NCNTRL_BITS_R-1:0];
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_NOTIFY) begin
			send_notify <= 1;
			send_notify_data <= mp_creg[31:0];
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_RESET) begin
			motion_reset <= 1;
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_WAIT_IDLE) begin
			if (m_cnt == 0) begin
				mp_cmd_end <= 1;
				len_fifo_rd_en <= 1;
			end
		end else if (mp_cmd == M_CMD_SET_EVENT) begin
			m_wait_mask <= mp_creg[NENDSTOP-1:0];
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_ENDSTOP_MASK) begin
			endstop_mask <= mp_creg[NENDSTOP-1:0];
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_ENDSTOP_POL) begin
			endstop_pol <= mp_creg[NENDSTOP-1:0];
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else if (mp_cmd == M_CMD_STOP) begin
			request_stop <= 1;
			mp_cmd_end <= 1;
			len_fifo_rd_en <= 1;
		end else begin
			/*
			 * unknwon command
			 */
			m_state_error <= 1;
		end
	end else if (mp_cmd_end) begin
		/*
		 * end of command stage. delayed by one clock so the recv
		 * len fifo can advance
		 */
		mp_cmd_end <= 0;
		mp_in_cmd <= 0;
	end

	/* endstop triggered, stop everything */
	if (endstop_int & endstop_mask) begin
		motion_reset <= 0;
		m_state_error <= 1;
	end

	/*
	 * motion control
	 */
	for (i = 0; i < NCNTRL; i = i + 1) begin
		if (do_step[i]) begin
			c_step[i + 1] <= ~c_step[i + 1];
		end
		do_step[i] <= 0;
	end
	if (m_cnt != 0) begin
		for (i = 0; i < NCNTRL; i = i + 1) begin
			m_snap[i] <= m_snap[i] + m_crackle[i];
			m_jerk[i] <= m_jerk[i] + m_snap[i];
			m_accel[i] <= m_accel[i] + m_jerk[i];
			m_velocity[i] <= m_velocity[i] + m_accel[i];
			/* do_step is the overflow of the addition */
			{ do_step[i], m_pos[i] } <=
				({ m_pos[i][REGBITS-1], m_pos[i] } +
			         { m_velocity[i][REGBITS-1], m_velocity[i] }) ^
				{ m_pos[i][REGBITS-1], { REGBITS { 1'b0 }}};
		end
	end
	if (motion_reset) begin
		for (i = 0; i < NCNTRL; i = i + 1) begin
			m_preload_jerk[i] <= 0;
			m_preload_snap[i] <= 0;
			m_preload_crackle[i] <= 0;
			m_jerk[i] <= 0;
			m_accel[i] <= 0;
			m_velocity[i] <= 0;
			do_step[i] <= 0;
			m_pos[i] <= 0;
			c_step[i + 1] <= 0;
		end
		for (i = 0; i < NSTEPDIR; i = i + 1) begin
			stepdir_routing[i] <= 0;
		end
		motion_reset <= 0;
	end
end

reg [2:0] notify_state = 0;
localparam N_IDLE = 0;
localparam N_WRITE_DATA1 = 1;
localparam N_WRITE_DATA2 = 2;
localparam N_WRITE_DATA3 = 3;
localparam N_WRITE_DATA4 = 4;
localparam N_WRITE_LEN = 5;
always @(posedge clk) begin
	if (send_ring_wr_en)
		send_ring_wr_en <= 0;
	if (send_fifo_wr_en)
		send_fifo_wr_en <= 0;
	if (send_notify) begin
		notify_state <= N_WRITE_DATA1;
	end else if (notify_state == N_WRITE_DATA1 && !send_ring_full) begin
		send_ring_data <= send_notify_data[31:24];
		send_ring_wr_en <= 1;
		notify_state <= N_WRITE_DATA2;
	end else if (notify_state == N_WRITE_DATA2 && !send_ring_full) begin
		send_ring_data <= send_notify_data[23:16];
		send_ring_wr_en <= 1;
		notify_state <= N_WRITE_DATA3;
	end else if (notify_state == N_WRITE_DATA3 && !send_ring_full) begin
		send_ring_data <= send_notify_data[15:8];
		send_ring_wr_en <= 1;
		notify_state <= N_WRITE_DATA4;
	end else if (notify_state == N_WRITE_DATA4 && !send_ring_full) begin
		send_ring_data <= send_notify_data[7:0];
		send_ring_wr_en <= 1;
		notify_state <= N_WRITE_LEN;
	end else if (notify_state == N_WRITE_LEN && !send_fifo_full) begin
		send_fifo_data <= 4;
		send_fifo_wr_en <= 1;
		notify_state <= N_IDLE;
	end
end
		
genvar gi;
generate
	for (gi = 0; gi < NCNTRL; gi = gi + 1) begin
		assign c_dir[gi + 1] = m_velocity[gi][REGBITS-1];
	end
endgenerate 

generate
for (gi = 0; gi < NSTEPDIR; gi = gi + 1) begin
	always @(*) begin
		step[gi] = c_step[stepdir_routing[gi]];
		dir[gi] = c_dir[stepdir_routing[gi]];
	end
end
endgenerate

assign debug = { mp_cmd, mp_in_cmd, running, m_state_error, 21'b0 };

endmodule
