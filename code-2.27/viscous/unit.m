% this is the bessel function for 1-j0(alpha*r)/j0(alpha*R0)
function out2=unit(r,R0,freq)
	omega=freq*2*pi;
	% check units for kinematic viscosity. ok.
	nu=1.0;
	alpha=sqrt(i*omega/nu);
	ar=alpha*r;
	aR0=alpha*R0;
	out2=1-(bessela(0,ar)/bessela(0,aR0));


