`timescale 1ns / 1ps
`default_nettype none

module pfs #(
	parameter BAUD = 9600,
	parameter NSTEPDIR = 6,
	parameter NCNTRL = 4
) (
	input clk_50mhz,

	output led,

	input rx1,
	output tx1,
	output cts1,

	input rx2,
	output tx2,
	output cts2,

	// stepper interface
	output en,
	output sck,
	output cs123,
	output cs456,
	output sdi,
	input sdo,
	output [NSTEPDIR-1:0] dir,
	output [NSTEPDIR-1:0] step,

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
	.send_ring_full(send_ring_full_1)
);

wire [31:0] motion_debug;
motion #(
	.LEN_BITS(LEN_BITS),
	.RECV_BUF_BITS(RECV_BUF_BITS),
	.NSTEPDIR(NSTEPDIR),
	.NCNTRL(NCNTRL)
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
control #(
	.LEN_BITS(LEN_BITS),
	.RECV_BUF_BITS(RECV_BUF_BITS),
	.NGPOUT(1),
	.NGPIN(1)
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

	.gpout({ en }),
	.gpin( { 1'b0 }),
	.sck(sck),
	.cs( { cs123, cs456 } ),
	.sdi(sdi),
	.sdo(sdo),

	/* debug */
	.debug(control_debug)
);

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
assign leds[223:184] = 0;

led7219 led7219_u(
	.clk(clk),
	.data(leds),
	.leds_out(leds_out),
	.leds_clk(leds_clk),
	.leds_cs(leds_cs)
);

endmodule
