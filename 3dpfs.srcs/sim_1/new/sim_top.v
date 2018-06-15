`timescale 1ns / 1ps

module sim_top();

reg clk_50mhz = 0;

wire leds_out;
wire leds_cs;
wire leds_clk;

main u_main(
	.clk_50mhz(clk_50mhz),

	.led(),

	.rx(),
	.tx(),

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

initial begin: B_clk
	integer i;
	for (i = 1; i < 100000; i = i + 1) begin
		clk_50mhz = 1;
		#10;
		clk_50mhz = 0;
		#10;
	end
end

endmodule
