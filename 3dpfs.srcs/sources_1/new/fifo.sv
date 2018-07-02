`timescale 1ns / 1ps
`default_nettype none

module fifo #(
	parameter DATA_WIDTH = 72,
	parameter ADDR_WIDTH = 6
) (
	input wire clk,
	input wire clr,
	// write side
	input wire [DATA_WIDTH - 1 : 0] din,
	input wire wr_en,
	output wire full,
	// read side
	output wire [DATA_WIDTH - 1 : 0] dout,
	input wire rd_en,
	output wire empty,

	// status
	output wire [ADDR_WIDTH - 1 : 0] elemcnt
);

localparam ADDRS = 1 << ADDR_WIDTH;
reg [DATA_WIDTH - 1 : 0] ram[ADDRS - 1 : 0];

reg [ADDR_WIDTH - 1 : 0] rdptr = 0;
reg [ADDR_WIDTH - 1 : 0] wrptr = 0;

wire [ADDR_WIDTH - 1 : 0] next_rdptr = rdptr + 1;
wire [ADDR_WIDTH - 1 : 0] next_wrptr = wrptr + 1;

assign empty = wrptr == rdptr;
assign full = next_wrptr == rdptr;
assign dout = ram[rdptr];
assign elemcnt = wrptr - rdptr;

always @(posedge clk) begin
	if (clr) begin
		rdptr <= 0;
		wrptr <= 0;
	end else begin
		if (rd_en && !empty) begin
			rdptr <= next_rdptr;
		end
		if (wr_en && !full) begin
			ram[wrptr] <= din;
			wrptr <= next_wrptr;
		end
	end
end

endmodule
