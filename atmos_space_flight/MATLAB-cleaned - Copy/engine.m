% Copyright Ashish Tewari (c) 2006
function [T, cT] = engine(alt, mach)
	% Mach number
	M=[0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25]; 
	
	% std. altitude (m)
	h=[0, 10, 20, 30, 36, 40, 50]*1000/3.28; 
	
	% Thrust (N)
	Thrust = [30, 21.5, 15, 10, 8, 6.5, 4; ...
		29, 21, 14.5, 9.8, 7.5, 6, 3.8; ...
		32, 22.5, 16, 10.5, 8.5, 7, 4.5; ...
		33, 28, 19, 12.5, 10, 8, 5; ...
		35, 29, 23.5, 16, 12.5, 10, 6; ...
		37, 31, 25.5, 21, 16, 13, 8.5; ...
		42.5, 35, 28, 22.5, 19.5, 15.5, 9.2; ...
		43.5, 38, 33, 25, 21.5, 17.5, 10.5; ...
		46, 39, 34, 28, 24.5, 19, 11.5; ...	
		48, 42, 35, 29, 26, 21.5, 13];
	Thrust = Thrust*1000*9.8/2.2;

	% (per hour)
	TSFC=[1.64, 1.66, 1.68, 1.7, 1.71, 1.71, 1.71; ...
		1.74, 1.76, 1.77, 1.78, 1.79, 1.79, 1.79; ...
		1.78, 1.79, 1.8, 1.815, 1.82, 1.82, 1.82; ...
		1.86, 1.8, 1.81, 1.82, 1.825, 1.825, 1.825; ...
		1.93, 1.84, 1.78, 1.79, 1.79, 1.79, 1.79; ...
		2, 1.9, 1.825, 1.76, 1.75, 1.75, 1.75; ...
		2.04, 1.96, 1.87, 1.79, 1.74, 1.74, 1.74; ...
		2.16, 2.05, 1.92, 1.84, 1.79, 1.79, 1.79; ...
		2.32, 2.14, 1.98, 1.88, 1.83, 1.83, 1.83; ...
		2.44, 2.26, 2.1, 1.97, 1.88, 1.88, 1.88];

	[X, Y] = meshgrid(h, M);
	T = interp2(X, Y, Thrust, alt, mach);
	cT = interp2(X, Y, TSFC, alt, mach);
end
