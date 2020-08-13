% Copyright Ashish Tewari (c) 2006
function p = Nstage(vf,beta,epsilon,alpha)
	p = 0.1;
	tol = 1e-9;

	N = size(beta,1);
	f = vf;

	for k = 1:N
		f = f+beta(k)*log(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
	end
	
	while abs(f) > tol
		f = vf;
		fp = 0;
		for k = 1:N
			f = f+beta(k)*log(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
			fp = fp+alpha(k)*beta(k)/(epsilon(k)+alpha(k)*(1-epsilon(k))*p);
		end
		d = -f/fp;
		p = p+d;
	end
end
