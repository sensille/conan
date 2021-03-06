`timescale 1ns / 1ps
`default_nettype none

module sim_top();

reg clk_50mhz = 0;

wire leds_out;
wire leds_cs;
wire leds_clk;

wire rx;
reg tx = 1;
wire rx2;
reg tx2 = 1;

wire cs123;
reg sdo;
wire sdi;
wire sck;

localparam BITRATE = 400000;
//localparam BITRATE = 9600;
localparam BITLENGTH = 1000000000 / BITRATE;

reg [7:0] endstop = 0;

pfs #(
	.BAUD(BITRATE),
	.NENDSTOP(8)
) u_pfs(
	.clk_50mhz(clk_50mhz),

	.led(),

	.rx1(tx),
	.tx1(rx),
	.cts1(),

	.rx2(tx2),
	.tx2(rx2),
	.cts2(),

	// stepper interface
	.en(),
	.sck(sck),
	.cs123(cs123),
	.cs456(),
	.sdi(sdi),
	.sdo(sdo),
	.dir(),
	.step(),

	.endstop(endstop),
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

task send2;
input [7:0] data;
begin : sendtask2
	integer i;

	tx2 = 0;	/* start bit */
	for (i = 0; i < 8; i = i + 1) begin
		#BITLENGTH tx2 = data[0];
		data = { 1'b0, data[7:1] };
	end
	#BITLENGTH tx2 = 1;	/* stop bit and idle */
	#BITLENGTH;
	#BITLENGTH;
end
endtask

initial begin: B_clk
	integer i;
	for (i = 0; i < 10000000; i = i + 1) begin
		clk_50mhz = 1;
		#10;
		clk_50mhz = 0;
		#10;
	end
end

initial begin: B_serial_data
	tx = 1;
	tx2 = 1;
	#1000;
	/* test gpout */
	send2(8'h80); send2(8'h70); send2(8'h02); send2(8'he9); send2(8'ha5);
	send2(8'h7e);
	send2(8'h80); send2(8'h70); send2(8'h03); send2(8'h29); send2(8'h64);
	send2(8'h7e);
	send2(8'h80); send2(8'h71); send2(8'h03); send2(8'hb9); send2(8'h65);
	send2(8'h7e);
	send2(8'h80); send2(8'h71); send2(8'h02); send2(8'h79); send2(8'ha4);
	send2(8'h7e);
	/* test SPI */
	send2(8'h01); send2(8'h80);
	send2(8'h01); send2(8'h02); send2(8'h03); send2(8'h04); send2(8'h05);
	send2(8'h11); send2(8'h12); send2(8'h13); send2(8'h14); send2(8'h15);
	send2(8'h21); send2(8'h22); send2(8'h23); send2(8'h24); send2(8'h25);
	send2(8'h71); send2(8'h19); send2(8'h7e);
	/* test motion controller */
	send(8'h80);
	send(8'ha0);
	send(8'h01);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h80); send(8'h60); send(8'h00); send(8'h01);
	send(8'h1e); send(8'he8); send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h80); send(8'h60); send(8'h02); send(8'h01);
	send(8'h7d); send(8'h5e); send(8'he9); send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h80);
	send(8'ha0);
	send(8'h01);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h01);
	send(8'h80);	/* LOAD reg 0 */
	send(8'h55);
	send(8'haa);
	send(8'h00);
	send(8'hb0);
	send(8'h7a);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h02);
	send(8'h71);
	send(8'h01);
	send(8'h00);
	send(8'h00);
	send(8'h3c);
	send(8'h33);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h03);
	send(8'h80);
	send(8'haa);
	send(8'h55);
	send(8'h01);
	send(8'h70);
	send(8'hb3);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h04);
	send(8'h71);
	send(8'h03);
	send(8'h00);
	send(8'h00);
	send(8'hfc);
	send(8'h1a);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h05);
	send(8'h80);	/* LOAD reg 0 */
	send(8'h55);
	send(8'haa);
	send(8'h00);
	send(8'h70);
	send(8'h8b);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h06);
	send(8'h71);
	send(8'h01);
	send(8'h00);
	send(8'h00);
	send(8'hfc);
	send(8'hc2);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h07);
	send(8'h61);
	send(8'h55);
	send(8'haa);
	send(8'h85);
	send(8'hef);
	send(8'h7e);
	#(BITLENGTH * 10);
	send(8'h08);
	send(8'h69);
	send(8'hff);
	send(8'hd2);
	send(8'hef);
	send(8'h7e);
	#(BITLENGTH * 10);

	/* start motion control */
	send2(8'h02); send2(8'h60); send2(8'h48); send2(8'h01); send2(8'h7e);
end

reg [119:0] data_out = 120'h391100006839110000683911000068;
always @(negedge sck) begin
	if (cs123 == 0) begin
		sdo <= data_out[119];
		data_out <= { data_out[118:0], 1'b0 };
	end
end
reg [119:0] data_in = 120'h0;
always @(posedge sck) begin
	if (cs123 == 0) begin
		data_in <= { data_in[118:0], sdi };
	end
end

initial begin: B_endstop
	integer i;
	for (i = 0; i < 10000; i = i + 1) begin
		endstop[0] = 1;
		#100;
		endstop[0] = 0;
		#100000;
	end
end

endmodule
