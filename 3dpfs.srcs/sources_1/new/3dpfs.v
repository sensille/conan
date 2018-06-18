`timescale 1ns / 1ps

module pfs(
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
wire [7:0] rx_data;
wire rx_ready;
reg [7:0] tx_data = 0;
reg tx_en = 0;
wire tx_is_transmitting;

parameter CLOCK_DIVIDE = 50; // clock rate (20Mhz) / (baud rate (100000) * 4)
uart #(
        .CLOCK_DIVIDE(CLOCK_DIVIDE)
) uart_u (
	.clk(clk),
	.rst(1'b0),
	.rx(rx_sync),
	.tx(tx),
	.transmit(tx_en),
	.tx_byte(tx_data),
	.received(rx_ready),
	.rx_byte(rx_data),
	.is_receiving(),
	.is_transmitting(tx_is_transmitting),
	.recv_error()
);

reg [7:0] recv_ring [RECV_BUF_SIZE-1:0];
reg [RECV_BUF_BITS-1:0] recv_wptr = 0;
reg [RECV_BUF_BITS-1:0] recv_wptr_fallback = 0;
reg [RECV_BUF_BITS-1:0] recv_rptr = 0;
reg recv_in_escape = 0;
reg [7:0] recv_seq = 0;

localparam RECV_ESC_CHAR = 8'h7d;
localparam RECV_EOF_CHAR = 8'h7e;

reg [7:0] recv_prev1 = 0;
reg [7:0] recv_prev2 = 0;
reg [7:0] crc16_in = 0;
reg [3:0] crc16_cnt = 0;
reg [15:0] crc16 = 0;
reg [7:0] recv_len = 0;
reg recv_len_pre1 = 0;
reg recv_len_pre2 = 0;
reg recv_len_pre3 = 0;
reg recv_error = 0;
reg recv_error_bad_state = 0;
reg recv_error_bad_crc = 0;
reg recv_error_bad_seq = 0;
reg recv_ack_sent = 0;
reg recv_send_ack = 0;
reg send_ack = 0;	/* signal for the send state machine */
reg [7:0] recv_last_seq = 0;

reg recv_len_wr_en = 0;

fifo #(
	.DATA_WIDTH(8)
) recv_len_fifo (
	.clk(clk),
	.clr(1'b0),

	// write side
	.din(recv_len),
	.wr_en(recv_len_wr_en),
	.full(),

	// read side
	.dout(),
	.rd_en(),
	.empty(),

	// status
	.elemcnt()
);

/*
 * packet receive state machine
 */
always @(posedge clk) begin
	recv_len_wr_en <= 0;
	if (rx_ready) begin
		if (rx_data == RECV_ESC_CHAR) begin
			if (recv_in_escape) begin
				recv_error <= 1;
				recv_error_bad_state <= 1;
				recv_send_ack <= 1;
			end else begin
				recv_in_escape <= 1;
			end
		end else if (rx_data == RECV_EOF_CHAR) begin
			if (recv_error) begin
				recv_wptr <= recv_wptr_fallback;
			end else if (recv_in_escape) begin
				recv_wptr <= recv_wptr_fallback;
				recv_error_bad_state <= 1;
				recv_send_ack <= 1;
			end else if (crc16 != { recv_prev2, recv_prev1 }) begin
				recv_wptr <= recv_wptr_fallback;
				recv_error_bad_crc <= 1;
				recv_send_ack <= 1;
			end else if (recv_seq != recv_last_seq + 1) begin
				recv_wptr <= recv_wptr_fallback;
				recv_error_bad_seq <= 1;
				recv_send_ack <= 1;
			end else begin
				recv_len_wr_en <= 1;
				recv_send_ack <= 1;
			end
			recv_error <= 0;
			recv_error_bad_state <= 0;
			recv_error_bad_crc <= 0;
			recv_error_bad_seq <= 0;
			recv_in_escape <= 0;
			recv_len <= 0;
			recv_len_pre1 <= 0;
			recv_len_pre2 <= 0;
			recv_len_pre3 <= 0;
			crc16 <= 0;
			recv_ack_sent <= 0;
			recv_last_seq <= recv_seq;
		end else if (recv_error) begin
			/* do nothing, wait for EOF */
		end else begin
			if (recv_in_escape)
				recv_prev1 <= rx_data ^ 8'h20;
			else
				recv_prev1 <= rx_data;
			recv_in_escape <= 0;
			recv_prev2 <= recv_prev1;
			recv_len_pre3 <= recv_len_pre2;
			recv_len_pre2 <= recv_len_pre1;
			recv_len_pre1 <= 1;
			if (recv_len_pre2 && !recv_len_pre3) begin
				recv_seq <= recv_prev2;
			end
			if (recv_len_pre3) begin
				recv_wptr <= recv_wptr + 1;
				recv_ring[recv_wptr] <= recv_prev2;
				recv_len <= recv_len + 1;
			end
			if (recv_len_pre2) begin
				crc16_in <= recv_prev2;
				crc16_cnt <= 8;
			end
		end
	end
	if (crc16_cnt != 0) begin
		crc16 <= { 1'b0, crc16[15:1] };
		crc16[15] <= crc16_in[0] ^ crc16[0];
		crc16[13] <= crc16[14] ^ crc16_in[0] ^ crc16[0];
		crc16[0] <= crc16[1] ^ crc16_in[0] ^ crc16[0];
		crc16_cnt <= crc16_cnt - 1;
		crc16_in <= { 1'b0, crc16_in[7:1] };
	end
	/*
	 * debounce: send only one ack per received frame
	 */
	if (recv_send_ack && !recv_ack_sent) begin
		recv_ack_sent <= 1;
		recv_send_ack <= 0;
		send_ack <= 1;
	end
	/*
	 * hold send_ack for only 1T
	 */
	if (send_ack) begin
		send_ack <= 0;
	end
end

localparam SEND_IDLE = 0;
localparam SEND_ACK_1 = 1;
localparam SEND_ACK_1_WAIT = 2;
localparam SEND_ACK_2 = 3;
localparam SEND_ACK_2_WAIT = 4;
reg [2:0] send_state = SEND_IDLE;
/*
 * packet send state machine
 */
always @(posedge clk) begin
	if (send_ack) begin
		send_state <= SEND_ACK_1;
	end
	if (!tx_is_transmitting) begin
		case (send_state)
			SEND_IDLE: begin
				tx_data <= 0;
			end
			SEND_ACK_1: begin
				tx_data <= {
					5'h00,
					recv_error_bad_state,
					recv_error_bad_crc,
					recv_error_bad_seq
				};
				tx_en <= 1;
				send_state <= SEND_ACK_1_WAIT;
			end
			SEND_ACK_1_WAIT: begin
				/* wait 1T for transmitting to go high */
				send_state <= SEND_ACK_2;
			end
			SEND_ACK_2: begin
				tx_data <= recv_last_seq;
				tx_en <= 1;
				send_state <= SEND_ACK_2_WAIT;
			end
			SEND_ACK_2_WAIT: begin
				/* wait 1T for transmitting to go high */
				send_state <= SEND_IDLE;
			end
		endcase
	end
	if (tx_en)
		tx_en <= 0;
end


/*
assign tx_data = rx_data ^ 8'h23;
assign tx_en = rx_ready;

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
*/


/*
// TODO: check channel

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
