`timescale 1ns / 1ps
`default_nettype none

module framing #(
	parameter BAUD = 9600,
	parameter RING_BITS = 10,
	parameter LEN_BITS = 8,
	parameter LEN_FIFO_BITS = 7,
	parameter HZ = 20000000
) (
	input wire clk,

	input wire rx,
	output wire tx,
	output wire cts,

	/*
	 * receive side
	 */
	/* len fifo output */
	output wire recv_fifo_empty,
	output wire [LEN_BITS-1:0] recv_fifo_dout,
	input wire recv_fifo_rd_en,

	/* ring buffer output */
	output wire [7:0] ring_data,
	input wire [RING_BITS-1:0] recv_rptr,

	/*
	 * send side
	 */
	input wire send_fifo_wr_en,
	input wire [LEN_BITS-1:0] send_fifo_data,
	output wire send_fifo_full,

	/* ring buffer input */
	input wire [7:0] send_ring_data,
	input wire send_ring_wr_en,
	output wire send_ring_full
);

localparam RING_SIZE = 1 << RING_BITS;

/*
 * sync rx into our clk domain
 */
reg rx_s1 = 0;
reg rx_sync = 0;
always @(posedge clk) begin
	rx_s1 <= rx;
	rx_sync <= rx_s1;
end

/*
 * UART instantiation
 */
wire [7:0] rx_data;
wire rx_ready;
reg [7:0] tx_data = 0;
reg tx_en = 0;
wire tx_transmitting;

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
	.is_transmitting(tx_transmitting),
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
reg [6:0] send_ack_seq = 0;
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

/*
 * send side
 */
wire [LEN_BITS-1:0] send_fifo_rd_data;
reg send_fifo_rd_en = 0;
wire send_fifo_empty;
fifo #(
	.DATA_WIDTH(LEN_BITS),
	.ADDR_WIDTH(LEN_FIFO_BITS)
) send_len_fifo (
	.clk(clk),
	.clr(1'b0),

	// write side
	.din(send_fifo_data),
	.wr_en(send_fifo_wr_en),
	.full(send_fifo_full),

	// read side
	.dout(send_fifo_rd_data),
	.rd_en(send_fifo_rd_en),
	.empty(send_fifo_empty),

	// status
	.elemcnt()
);

/*
 * sending data from the fifo takes precedence, as each data packet also
 * carries an ack seq number.
 */
localparam SEND_IDLE = 0;
localparam SEND_DATA = 1;
localparam SEND_CRC1 = 2;
localparam SEND_CRC2 = 3;
localparam SEND_EOF = 4;
reg [2:0] send_state = 0;
reg do_send = 0;
reg send_ack_requested = 0;
reg send_in_escape = 0;
reg [3:0] send_crc16_cnt = 0;
reg [7:0] send_crc16_in = 0;
reg [15:0] send_crc16 = 0;
reg [7:0] send_byte = 0;
reg [7:0] send_ring [RING_SIZE-1:0];
reg [RING_BITS-1:0] send_rptr = 0;
reg [RING_BITS-1:0] send_wptr = 0;
reg [LEN_BITS-1:0] send_len = 0;
assign send_ring_full = (send_wptr + 1) == send_rptr;

always @(posedge clk) begin
	if (send_ring_wr_en && !send_ring_full) begin
		send_ring[send_wptr] <= send_ring_data;
		send_wptr <= send_wptr + 1;
	end
end

always @(posedge clk) begin
	if (send_ack)
		send_ack_requested <= 1;
	if (send_fifo_rd_en)
		send_fifo_rd_en <= 0;
	if (tx_en)
		tx_en <= 0;
	if (send_crc16_cnt != 0) begin
		send_crc16 <= { 1'b0, send_crc16[15:1] };
		send_crc16[15] <= send_crc16_in[0] ^ send_crc16[0];
		send_crc16[13] <= send_crc16[14] ^ send_crc16_in[0] ^
			send_crc16[0];
		send_crc16[0] <= send_crc16[1] ^ send_crc16_in[0] ^
			send_crc16[0];
		send_crc16_cnt <= send_crc16_cnt - 1;
		send_crc16_in <= { 1'b0, send_crc16_in[7:1] };
	end
	if (do_send && !tx_transmitting && !tx_en) begin
		if (send_in_escape) begin
			tx_data <= send_byte ^ 8'h20;
			do_send <= 0;
			send_in_escape <= 0;
		end else begin
			if (send_byte == RECV_ESC_CHAR ||
			    send_byte == RECV_EOF_CHAR) begin
				tx_data <= RECV_ESC_CHAR;
				send_in_escape <= 1;
			end else begin
				tx_data <= send_byte;
				do_send <= 0;
			end
		end
		tx_en <= 1;
	end else if (tx_transmitting || tx_en) begin
		/* do nothing, wait for send to complete */
	end else if (send_state == SEND_IDLE &&
	    (send_ack_requested || !send_fifo_empty)) begin
		if (send_fifo_empty)
			send_len <= 0;
		else
			send_len <= send_fifo_data;
		send_ack_requested <= 0;
		send_state <= SEND_DATA;
		send_crc16 <= 0;
		send_byte <= { (send_ack_error_state != 0) ? 1'b1 : 1'b0,
		               send_ack_seq[6:0] };
		send_crc16_in <= { (send_ack_error_state != 0) ? 1'b1 : 1'b0,
		               send_ack_seq[6:0] };
		send_crc16_cnt <= 8;
		do_send <= 1;
		send_fifo_rd_en <= 1;
	end else if (send_state == SEND_DATA) begin
		if (send_len != 0) begin
			send_byte <= send_ring[send_rptr];
			send_crc16_in <= send_ring[send_rptr];
			send_rptr <= send_rptr + 1;
			send_len <= send_len - 1;
			do_send <= 1;
			send_crc16_cnt <= 8;
		end else begin
			send_state <= SEND_CRC1;
		end
	end else if (send_state == SEND_CRC1) begin
		send_byte <= send_crc16[15:8];
		do_send <= 1;
		send_state <= SEND_CRC2;
	end else if (send_state == SEND_CRC2) begin
		send_byte <= send_crc16[7:0];
		do_send <= 1;
		send_state <= SEND_EOF;
	end else if (send_state == SEND_EOF) begin
		tx_data <= RECV_EOF_CHAR;
		tx_en <= 1;
		send_state <= SEND_IDLE;
	end
end

endmodule
