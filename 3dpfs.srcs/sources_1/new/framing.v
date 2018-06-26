`timescale 1ns / 1ps

module framing #(
	parameter BAUD = 9600,
	parameter RING_BITS = 10,
	parameter LEN_BITS = 8,
	parameter LEN_FIFO_BITS = 7,
	parameter HZ = 20000000
) (
	input clk,

	input rx,
	output tx,
	output cts,

	/* len fifo output */
	output recv_fifo_empty,
	output [LEN_BITS-1:0] recv_fifo_dout,
	output recv_fifo_rd_en,

	/* ring buffer output */
	output [7:0] ring_data,
	input [RING_BITS-1:0] recv_rptr
);

localparam RING_SIZE = 1 << RING_BITS;

/*
 * sync rx into our clk domain
 */
reg rx_s1 = 0;
reg rx_sync = 0;
always @(posedge clk) begin
	rx_s1 = rx;
	rx_sync = rx_s1;
end

/*
 * UART instantiation
 */
wire [7:0] rx_data;
wire rx_ready;
reg [7:0] tx_data = 0;
reg tx_en = 0;
wire tx_is_transmitting;

localparam CLOCK_DIVIDE = HZ / BAUD / 4; /* 4 phases per bit */
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

reg [7:0] recv_ring [RING_SIZE-1:0];
reg [RING_BITS-1:0] recv_wptr = 0;
reg [RING_BITS-1:0] recv_wptr_fallback = 0;
reg recv_in_escape = 0;
reg [7:0] recv_seq = 0;
assign ring_data = recv_ring[recv_rptr];

localparam RECV_ESC_CHAR = 8'h7d;
localparam RECV_EOF_CHAR = 8'h7e;

localparam RECV_ERROR_OK = 0;
localparam RECV_ERROR_BAD_STATE = 1;
localparam RECV_ERROR_BAD_CRC = 2;
localparam RECV_ERROR_BAD_SEQ = 3;

reg [1:0] recv_error_state = RECV_ERROR_OK;

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
reg recv_ack_sent = 0;
reg recv_send_ack = 0;
reg [6:0] recv_last_seq = 0;
reg recv_eof = 0;
reg [7:0] recv_len_wr_data = 0;

reg recv_len_wr_en = 0;

wire recv_len_fifo_full;

assign cts = (recv_rptr == recv_wptr + 1) || recv_len_fifo_full;

fifo #(
	.DATA_WIDTH(LEN_BITS),
	.ADDR_WIDTH(LEN_FIFO_BITS)
) recv_len_fifo (
	.clk(clk),
	.clr(1'b0),

	// write side
	.din(recv_len_wr_data),
	.wr_en(recv_len_wr_en),
	.full(recv_len_fifo_full),

	// read side
	.dout(recv_fifo_dout),
	.rd_en(recv_fifo_rd_en),
	.empty(recv_fifo_empty),

	// status
	.elemcnt()
);

reg send_ack = 0;	/* signal for the send state machine */
reg [1:0] send_ack_error_state = 0;
reg [7:0] send_ack_seq = 0;
/*
 * packet receive state machine
 */
always @(posedge clk) begin
	recv_len_wr_en <= 0;
	if (rx_ready) begin
		if (rx_data == RECV_ESC_CHAR) begin
			if (recv_in_escape) begin
				recv_error <= 1;
				recv_error_state <= RECV_ERROR_BAD_STATE;
				recv_send_ack <= 1;
			end else begin
				recv_in_escape <= 1;
			end
		end else if (rx_data == RECV_EOF_CHAR) begin
			if (recv_error) begin
				recv_wptr <= recv_wptr_fallback;
			end else if (recv_in_escape) begin
				recv_wptr <= recv_wptr_fallback;
				recv_error_state <= RECV_ERROR_BAD_STATE;
			end else if (crc16 != { recv_prev2, recv_prev1 }) begin
				recv_wptr <= recv_wptr_fallback;
				recv_error_state <= RECV_ERROR_BAD_CRC;
			end else if (recv_seq[7] == 0 &&
			             recv_seq[6:0] != recv_last_seq + 1) begin
				recv_wptr <= recv_wptr_fallback;
				recv_error_state <= RECV_ERROR_BAD_SEQ;
			end else begin
				if (recv_len != 0)
					recv_len_wr_en <= 1;
				recv_len_wr_data <= recv_len;
				recv_wptr_fallback <= recv_wptr;
			end
			recv_in_escape <= 0;
			recv_len <= 0;
			recv_len_pre1 <= 0;
			recv_len_pre2 <= 0;
			recv_len_pre3 <= 0;
			crc16 <= 0;
			recv_last_seq <= recv_seq[6:0];
			recv_eof <= 1;
			recv_send_ack <= 1;
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
		/* values passed to send state machine */
		send_ack <= 1;
		send_ack_error_state <= recv_error_state;
		send_ack_seq <= recv_last_seq;
	end

	/*
	 * hold send_ack for only 1T
	 */
	if (send_ack) begin
		send_ack <= 0;
	end

	/*
	 * reset error state on frame end
	 */
	if (recv_eof) begin
		recv_eof <= 0;
		recv_error <= 0;
		recv_error_state <= RECV_ERROR_OK;
		recv_ack_sent <= 0;
		recv_send_ack <= 0;
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
					6'h00,
					send_ack_error_state
				};
				tx_en <= 1;
				send_state <= SEND_ACK_1_WAIT;
			end
			SEND_ACK_1_WAIT: begin
				/* wait 1T for transmitting to go high */
				send_state <= SEND_ACK_2;
			end
			SEND_ACK_2: begin
				tx_data <= send_ack_seq;
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

endmodule
