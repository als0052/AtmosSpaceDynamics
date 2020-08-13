% TODO: figure out units
% TODO: break up the code in lines 47-63 into separate function, call for each orbit
function [r, v] = relative(orb1, orb2, t, mu)
	%RELATIVE calculates the relative position and velocity of two elliptical orbits
	%
	% [r, v] = relative(orb1, orb2, t)
	%
	% Inputs:
	%   orb1: The elliptical orbit of the target spacecraft
	%   orb2: The elliptical orbit of the object (?) spacecraft
	%   t: 
	%   mu: Standard gravitatonial parameter. Default is 398600.4 km^3/s^2
	%
	% Outputs:
	%   r: The relative position (UNITS?)
	%   v: The relative velocity (UNITS?)
	% 
	% Each orbit (orb1, orb2) is expected to be a struct with the following
	% attributes: (the classic orbital elements)
	%   a: semi-major axis (UNITS?)
	%   e: eccentricity (dimensionless)
	%   i: orbit inclination (radians?)
	%   w: omega; argument of periapsis (UNITS?)
	%   Om: Omega; right ascension of the ascending node (UNITS?)
	%   tau: time of periapsis (UNITS?)
	% 
	% Requirements:
	%   This program calls trajE.m
	% 
	% See also:
	%   Page 145
	% 
	% Copyright Ashish Tewari (c) 2006
	% Modified by Andrew Smelser
	% 08-13-2020
	
	if nargin < 2 || isempty(orb1) || isempty(orb2)
		error("Both orbits are required inputs\n")
	end
	if nagrin < 3 || isempty(t)
		error("<???>, t, is a required input\n")
	end
	if nargin < 4 || isempty(mu)
		mu = 398600.4;	% km^3/s^2
	end
	
	a = orb1(1);
	e = orb1(2);
	i = orb1(3);
	w = orb1(4);
	Om = orb1(5);
	tau = orb1(6);
	
	n = sqrt(mu/a^3);
	M = -n*tau;
	E = kepler(e, M);
	r0 = a*(1 - e*cos(E));
	R0 = a*[cos(E) - e; sqrt(1 - e^2)*sin(E); 0];
	V0 = sqrt(mu*a)*[-sin(E); sqrt(1 - e^2)*cos(E); 0]/r0;
	[Rt, Vt] = trajE(mu, 0, R0, V0, t);
	C = rotation(i, Om, w);
	Rt = C*Rt;
	Vt = C*Vt;
	
	a = orb2(1);
	e = orb2(2);
	i = orb2(3);
	w = orb2(4);
	Om = orb2(5);
	tau = orb2(6);
	
	n = sqrt(mu/a^3);
	M = -n*tau;
	E = kepler(e, M);
	r0 = a*(1 - e*cos(E));
	R0 = a*[cos(E) - e; sqrt(1 - e^2)*sin(E); 0];
	V0 = sqrt(mu*a)*[-sin(E); sqrt(1 - e^2)*cos(E); 0]/r0;
	[R, V] = trajE(mu, 0, R0, V0, t);
	C = rotation(i, Om, w);
	R = C*R;
	V = C*V;
	r = R - Rt;
	v = V - Vt;
	rt = norm(Rt);
	
	lat = asin(dot(Rt, [0; 0; 1])/rt);
	slon = dot(Rt, [0; 1; 0])/(rt*cos(lat));
	clon = dot(Rt, [1; 0; 0])/(rt*cos(lat));
	long = atan(slon/clon);

	if slon < 0 && clon > 0
		long = asin(slon);
	elseif slon > 0 && clon < 0
		long = acos(clon);
	end
	
	CLH = INtoLH(lat, long);
	r = CLH*r;
	v = CLH*v;
end
