% TODO: figure out units
function [R, V] = trajE(mu, t0, R0, V0, t, tol)
	%trajE Determine the elliptical trajectory
	%
	% [R, V] = trajE(mu, t0, R0, V0, t)
	%
	% Inputs:
	%   mu: Standard gravitational perameter (UNITS?)
	%   t0: Initial (starting) time (seconds?)
	%   R0: Initial position vector (UNITS?)
	%   V0: Initial velocity vector (UNITS?)
	%   t: Final (ending) time (seconds?)
	%   tol: Tolerance value. Default is 1e-10
	%
	% Outputs:
	%   R: Position vector (UNITS?)
	%   V: Velocity vector (UNITS?)
	%
	% See also:
	%   Page 107
	% 
	% Copyright Ashish Tewari (c) 2006
	% Modified by Andrew Smelser
	% 08-13-2020
	
	if nargin < 1 || isempty(mu)
		error("Standard gravitational parameter, mu, is a required input\n")
	end
	if nargin < 2 || isempty(t0)
		t0 = 0;
	end
	if nargin < 3 || isempty(R0)
		error("Initial position vector, R0, is a required input\n")
	elseif length(R0) ~= 3 || ~isvector(R0)
		error("Initial position vector, R0, must be a vector of length 3\n")
	end
	if nargin < 4 || isempty(V0)
		error("Initial velocity vector, V0, is a required input\n")
	elseif length(V0) ~= 3 || ~isvector(V0)
		error("Initial velocity vector, V0, must be a vector of length 3\n")
	end
	if nargin < 5 || isempty(t)
		error("Final time, t, is a required input\n")
	end
	if nargin < 6 || isempty(tol)
		tol = 1e-10;
	end
	
	r0 = norm(R0);
	v0 = norm(V0);
	alpha = dot(R0, V0);
	H = cross(R0, V0);
	h = norm(H);
	p = h^2/mu;
	ecv = cross(V0, H)/mu - R0/r0;
	e = norm(ecv);
	ecth0 = p/r0 - 1;
	esth0 = norm(cross(ecv, R0))/r0;

	if abs(ecth0) >= tol
		th0 = atan(esth0/ecth0);
		if ecth0 < 0
			if esth0 >= 0
				th0 = th0 + pi;
			end
		elseif esth0 < 0
			th0 = th0 + 2*pi;
		end
	elseif esth0 >= 0
		th0 = pi/2;
	else
		th0 = 3*pi/2;
	end

	ainv = -(v0^2)/mu + 2/r0;

	if abs(1-e) > tol
		a = 1/ainv;
		if e < 1
			n = sqrt(mu/a^3);
			E0 = 2*atan(sqrt((1 - e)/(1 + e))*tan(0.5*th0));
			tau = t0 + (-E0 + e*sin(E0))/n;
			M = n*(t - tau);
			E = kepler(e,M);
			r = a*(1 - e*cos(E));
		else
			n = sqrt(-mu/a^3);
		end
	else
		e = 1;
	end

	f = 1 + a*(cos(E - E0) - 1)/r0;
	g = a*alpha*(1 - cos(E-E0))/mu + r0*sqrt(a/mu)*sin(E - E0);
	fd = -sqrt(mu*a)*(sin(E - E0))/(r*r0);
	gd = 1 + a*(cos(E - E0) - 1)/r;
	R = f*R0 + g*V0;
	V = fd*R0 + gd*V0;
end
