
%   ok, here's the scoop
%   there are three other .m files that get used -
%   unit, corey and coreyfile.

%   unit is the bessel function for the heart of gamma, 
%       i.e.  1-(j0(a*r)/j0(a*R)).
%   corey is the discretized calculation of the integral of unit over 
%   r from 0 to R0.  it is as good as quadrature method for integrating 
%   for discretization sufficiently high (i.e. n > 50).
%   and it is also more versatile - you can give it 
%   corey(R0,freq,n) where n is the discretization.
%   and it will run unit for whatever channel height (R0) and input
%   frequency that you like.  therefore, to build tables in AKOUEIN
%   for the viscous term, gamma, we run corey at n > 50 for R0=.03
%   (that's the tunnel) and R0=.06 (the sulcus) for all the input
%   frequencies in the free world.  the output of tables is a string
%   of real and imaginary terms
%   a table for the sulcus and a table for the tunnel
%
%   there is nothing sacred about these tables and fi they are missing a
%   frequency value that we use later, just run corey again as
%      corey(.03,freq,75) 
%   the tables are easy to generate, obviously, using table.m, but 
%   compulsively making sure that every queero input paradigm gets
%   automatically generated is madness, i tell you, madness.
%   oh and coreyfile is just a call to corey with some file writing

function table(n)
	fpss=fopen('tabless.test','w');
	fptc=fopen('tabletc.test','w');
	fprintf(fptc,'tunnel of corti - R0=.03 mm \n ');
	fprintf(fptc,'frequency\treal\t\timag\t\tmagni\t\tphase \n ');
	fprintf(fpss,'spiral sulcus - R0=.06 mm \n ');
	fprintf(fpss,'frequency\treal\t\timag\t\tmagni\t\tphase \n ');
	discretization=n
	freq=300;
	coreyfile(fptc,.03,freq,discretization)
	coreyfile(fpss,.06,freq,discretization)
	% get gamma for frequencies 1000 to 12,000 
	for i=1:24
		coreyfile(fptc,.03,i*500,discretization)
		coreyfile(fpss,.06,i*500,discretization)
	end;
	fclose(fpss);
	fclose(fptc);



