`timescale 1ns / 1ps
`default_nettype none

module pfs #(
	parameter BAUD = 115200,
	parameter NSTEPDIR = 6,
	parameter NCNTRL = 4,
	parameter NENDSTOP = 8
) (
	input wire clk_50mhz,

	output wire led,

	input wire rx1,
	output wire tx1,
	output wire cts1,

	input wire rx2,
	output wire tx2,
	output wire cts2,

	// stepper interface
	output wire en,
	output wire sck,
	output wire cs123,
	output wire cs456,
	output wire sdi,
	input wire sdo,
	output wire [NSTEPDIR-1:0] dir,
	output wire [NSTEPDIR-1:0] step,

	// endstops
	input wire [NENDSTOP-1:0] endstop,

	// debug
	output wire debug1,
	output wire debug2,
	output wire debug3,
	output wire debug4,

	//
	output wire leds_out,
	output wire leds_clk,
	output wire leds_cs
);

localparam LEN_BITS = 8;
localparam RECV_BUF_BITS = 10;
localparam RECV_BUF_SIZE = (1 << 10);
localparam HZ = 20000000;

assign debug1 = rx1;
assign debug2 = tx1;
assign debug3 = rx2;
assign debug4 = tx2;

wire clk;
clk u_clk(
	.clk_50mhz(clk_50mhz),
	.clk(clk)		// 20 MHz
);

reg [31:0] led_cnt = 0;
always @(posedge clk)
	led_cnt = led_cnt + 1;
assign led = led_cnt[22];

/*
 * UART 1 instantiation
 */
wire len_fifo_empty_1;
wire [7:0] len_fifo_data_1;
wire len_fifo_rd_en_1;
wire [7:0] ring_data_1;
wire [9:0] recv_rptr_1;
wire send_fifo_full_1;
wire send_fifo_wr_en_1;
wire [LEN_BITS-1:0] send_fifo_data_1;
wire send_ring_full_1;
wire send_ring_wr_en_1;
wire [7:0] send_ring_data_1;
framing #(
	.BAUD(BAUD),
	.RING_BITS(10),
	.LEN_BITS(LEN_BITS),
	.LEN_FIFO_BITS(7),
	.HZ(HZ)
) u_framing_1 (
	.clk(clk),

	.rx(rx1),
	.tx(tx1),
	.cts(cts1),

	/* len fifo output */
	.recv_fifo_empty(len_fifo_empty_1),
	.recv_fifo_dout(len_fifo_data_1),
	.recv_fifo_rd_en(len_fifo_rd_en_1),

	/* ring buffer output */
	.ring_data(ring_data_1),
	.recv_rptr(recv_rptr_1),

	/* send len fifo */
	.send_fifo_wr_en(send_fifo_wr_en_1),
	.send_fifo_data(send_fifo_data_1),
	.send_fifo_full(send_fifo_full_1),

	/* send ring */
	.send_ring_data(send_ring_data_1),
	.send_ring_wr_en(send_ring_wr_en_1),
	.send_ring_full(send_ring_full_1),

	.clr(clear_error)
);

wire [31:0] motion_debug;
wire running;
wire clear_error;
wire request_stop;
motion #(
	.LEN_BITS(LEN_BITS),
	.RECV_BUF_BITS(RECV_BUF_BITS),
	.NSTEPDIR(NSTEPDIR),
	.NCNTRL(NCNTRL),
	.NENDSTOP(NENDSTOP)
) u_motion (
	.clk(clk),

	/* len fifo input */
	.len_fifo_empty(len_fifo_empty_1),
	.len_fifo_data(len_fifo_data_1),
	.len_fifo_rd_en(len_fifo_rd_en_1),

	/* ring buffer input */
	.recv_data(ring_data_1),
	.recv_rptr(recv_rptr_1),

	/* send len fifo */
	.send_fifo_wr_en(send_fifo_wr_en_1),
	.send_fifo_data(send_fifo_data_1),
	.send_fifo_full(send_fifo_full_1),

	/* send ring */
	.send_ring_data(send_ring_data_1),
	.send_ring_wr_en(send_ring_wr_en_1),
	.send_ring_full(send_ring_full_1),

	.dir(dir),
	.step(step),

	.running(running),
	.clear_error(clear_error),
	.request_stop(request_stop),
	.endstop(endstop),

	/* debug */
	.debug(motion_debug)
);

/*
 * UART 2 instantiation
 */
wire len_fifo_empty_2;
wire [7:0] len_fifo_data_2;
wire len_fifo_rd_en_2;
wire [7:0] ring_data_2;
wire [9:0] recv_rptr_2;

wire send_fifo_full_2;
wire send_fifo_wr_en_2;
wire [LEN_BITS-1:0] send_fifo_data_2;
wire send_ring_full_2;
wire send_ring_wr_en_2;
wire [7:0] send_ring_data_2;
framing #(
	.BAUD(BAUD),
	.RING_BITS(10),
	.LEN_BITS(LEN_BITS),
	.LEN_FIFO_BITS(7),
	.HZ(HZ)
) u_framing_2 (
	.clk(clk),

	.rx(rx2),
	.tx(tx2),
	.cts(cts2),

	/* len fifo output */
	.recv_fifo_empty(len_fifo_empty_2),
	.recv_fifo_dout(len_fifo_data_2),
	.recv_fifo_rd_en(len_fifo_rd_en_2),

	/* ring buffer output */
	.ring_data(ring_data_2),
	.recv_rptr(recv_rptr_2),

	/* send len fifo */
	.send_fifo_wr_en(send_fifo_wr_en_2),
	.send_fifo_data(send_fifo_data_2),
	.send_fifo_full(send_fifo_full_2),

	/* send ring */
	.send_ring_data(send_ring_data_2),
	.send_ring_wr_en(send_ring_wr_en_2),
	.send_ring_full(send_ring_full_2)
);

wire [31:0] control_debug;
localparam NGPOUT = 8;
wire [NGPOUT-1:0] gpout;
control #(
	.LEN_BITS(LEN_BITS),
	.RECV_BUF_BITS(RECV_BUF_BITS),
	.NGPOUT(NGPOUT),
	.NGPIN(NENDSTOP)
) u_control (
	.clk(clk),

	/* len fifo input */
	.len_fifo_empty(len_fifo_empty_2),
	.len_fifo_data(len_fifo_data_2),
	.len_fifo_rd_en(len_fifo_rd_en_2),

	/* ring buffer input */
	.recv_data(ring_data_2),
	.recv_rptr(recv_rptr_2),

	/* send len fifo */
	.send_fifo_wr_en(send_fifo_wr_en_2),
	.send_fifo_data(send_fifo_data_2),
	.send_fifo_full(send_fifo_full_2),

	/* send ring */
	.send_ring_data(send_ring_data_2),
	.send_ring_wr_en(send_ring_wr_en_2),
	.send_ring_full(send_ring_full_2),

	.gpout(gpout),
	.gpin(endstop),
	.sck(sck),
	.cs( { cs123, cs456 } ),
	.sdi(sdi),
	.sdo(sdo),

	.running(running),
	.request_stop(request_stop),
	.clear_error(clear_error),

	/* debug */
	.debug(control_debug)
);

assign en = gpout[0];

wire [255:0] leds;
assign leds[47:0] = { len_fifo_empty_1, len_fifo_data_1, len_fifo_rd_en_1,
	ring_data_1, recv_rptr_1, send_fifo_full_1, send_fifo_wr_en_1,
	send_fifo_data_1, send_ring_full_1, send_ring_wr_en_1,
	send_ring_data_1 };
assign leds[95:48] = { len_fifo_empty_2, len_fifo_data_2, len_fifo_rd_en_2,
	ring_data_2, recv_rptr_2, send_fifo_full_2, send_fifo_wr_en_2,
	send_fifo_data_2, send_ring_full_2, send_ring_wr_en_2,
	send_ring_data_2 };
assign leds[127:96] = motion_debug;
assign leds[159:128] = control_debug;
assign leds[255:224] = led_cnt;

/* 24 signals */
assign leds[183:160] = { rx1, tx1, cts1, rx2, tx2, cts2, en,
	sck, cs123, cs456, sdi, sdo, dir, step };
assign leds[223:216] = gpout;
assign leds[191:184] = endstop;
assign leds[215:192] = 0;

led7219 led7219_u(
	.clk(clk),
	.data(leds),
	.leds_out(leds_out),
	.leds_clk(leds_clk),
	.leds_cs(leds_cs)
);

endmodule
