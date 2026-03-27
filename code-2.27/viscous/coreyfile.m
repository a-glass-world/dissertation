
function coreyfile(fp,R0,freq,n)
	out=corey(R0,freq,n)
	fprintf(fp,'%e ',freq);
	fprintf(fp,'%s\t ',out);
	fprintf(fp,'\n ');



