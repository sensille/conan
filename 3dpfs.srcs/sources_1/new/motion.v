`timescale 1ns / 1ps

module motion #(
	parameter LEN_BITS = 8,
	parameter RECV_BUF_BITS = 10
) (
	input clk,

	/* len fifo input */
	input len_fifo_empty,
	input [LEN_BITS-1:0] len_fifo_data,
	output reg len_fifo_rd_en,

	/* ring buffer input */
	input [7:0] recv_data,	/* data at rptr */
	output reg [RECV_BUF_BITS-1:0] recv_rptr = 0,

	/* debug output */
	output [31:0] debug
);

reg [31:0] r2nd_cnt = 0;
reg [31:0] r2nd_cnt_preload = 0;
reg r2nd_in_read = 0;
reg [7:0] r2nd_read_len = 0;
always @(posedge clk) begin
	len_fifo_rd_en <= 0;
	if (r2nd_cnt != 0) begin
		r2nd_cnt <= r2nd_cnt - 1;
	end else if (!len_fifo_empty && r2nd_in_read == 0) begin
		r2nd_cnt_preload <= 0;
		r2nd_read_len <= len_fifo_data;
		r2nd_in_read <= 1;
	end else if (r2nd_read_len != 0) begin
		r2nd_read_len <= r2nd_read_len - 1;
		r2nd_cnt_preload <= { r2nd_cnt_preload[23:0], recv_data };
		recv_rptr <= recv_rptr + 1;
	end else if (r2nd_in_read) begin
		r2nd_in_read <= 0;
		len_fifo_rd_en <= 1;
		r2nd_cnt <= r2nd_cnt_preload;
	end else begin
		len_fifo_rd_en <= 0;
	end
end

assign debug = r2nd_cnt;

endmodule
