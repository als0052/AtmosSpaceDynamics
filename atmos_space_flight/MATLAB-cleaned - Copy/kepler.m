% TODO: figure out units
function [E, iter] = kepler(e, M, tol, maxIter, verbose)
	%kepler Solving for Kepler's equation
	%
	% E = kepler(e, M)
	%
	% Inputs:
	%   e: Eccentricity (dimensionless)
	%   M: Mean anomaly (UNITS?)
	%   tol: Tolerance value. Default is 1e-10
	%   maxIter: The maximum number of allowable iterations. Default is 10.
	%   verbose: Boolean which determines whether or not to print the results to
	%            the console. Defaut is False
	%
	% Outputs:
	%   E: Eccentric anomaly (UNITS?)
	%   iter: The number of iterations it took for the function to converge
	%
	% See also:
	%   Page 105
	%
	% Copyright Ashish Tewari (c) 2006
	% Modified by Andrew Smelser
	% 08-13-2020
	
	if nargin < 2 || isempty(e) || isempty(M)
		error("Eccentricity, e, and Mean anomaly, M, are required inputs\n");
	end
	if nargin < 3 || isempty(tol)
		tol = 1e-10;
	end
	if nargin < 4 || isempty(maxIter)
		maxIter = 10;
	end
	if nargin < 5 || isempty(verbose)
		verbose = false;
	end
	
	iter = 0;
	B = cos(e) - (pi/2 - e)*sin(e);
	E = M + e*sin(M)/(B + M*sin(e));
	fE = E - e*sin(E) - M;
	fpE = 1 - e*cos(E);
	dE = -fE/fpE;
	E = E + dE;

	while abs(fE) >= tol
		fE = E - e*sin(E) - M;
		fpE = 1 - e*cos(E);
		dE = -fE/fpE;
		E = E + dE;
		iter = iter + 1;
		if iter > maxIter
			fprintf("kepler.m\n\tMaximum allowable iterations exceeded\n")
			break
		end
	end
	
	if verbose
		fprintf("kepler.m\n");
		fprintf("E:%11.4f\n", E)
		fprintf("Iterations:%2d\n\n", iter)
	end
end
