`timescale 1ns / 1ps

module sim_top();

reg clk_50mhz = 0;

wire leds_out;
wire leds_cs;
wire leds_clk;

wire rx;
reg tx = 1;

localparam BITRATE = 100000;
localparam BITLENGTH = 1000000000 / BITRATE;

pfs #(
	.BAUD(BITRATE)
) u_pfs(
	.clk_50mhz(clk_50mhz),

	.led(),

	.rx(tx),
	.tx(rx),
	.cts(),

	// stepper interface
	.en(),
	.sck(),
	.cs123(),
	.cs456(),
	.sdi(),
	.sdo(),
	.dir(),
	.step(),

	// debug
	.debug1(),
	.debug2(),
	.debug3(),
	.debug4(),

	//
	.leds_out(leds_out),
	.leds_clk(leds_clk),
	.leds_cs(leds_cs)
);

task send;
input [7:0] data;
begin : sendtask
	integer i;

	tx = 0;	/* start bit */
	for (i = 0; i < 8; i = i + 1) begin
		#BITLENGTH tx = data[0];
		data = { 1'b0, data[7:1] };
	end
	#BITLENGTH tx = 1;	/* stop bit and idle */
	#BITLENGTH;
	#BITLENGTH;
end
endtask

initial begin: B_clk
	integer i;
	for (i = 0; i < 1000000; i = i + 1) begin
		clk_50mhz = 1;
		#10;
		clk_50mhz = 0;
		#10;
	end
end

initial begin: B_serial_data
	tx = 1;
	#1000;
	send(8'h3a);
	send(8'h5b);
	send(8'h70);
	send(8'h7d);
	send(8'h5d);
	send(8'h7d);
	send(8'h5e);
	send(8'hd9);
	send(8'h75);
	send(8'h0a);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h80);
	send(8'ha0);
	send(8'h01);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h01);
	send(8'h55);
	send(8'h02);
	send(8'h51);
	send(8'hef);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h02);
	send(8'haa);
	send(8'h01);
	send(8'ha0);
	send(8'h1e);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h7d);
	send(8'h7d);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h02);
	send(8'hc0);
	send(8'hc2);
	send(8'h7e);
end

endmodule
