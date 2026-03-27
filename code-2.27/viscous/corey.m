function out=corey(R0,freq,n)

% calculate   1 - j0(alpha*r)/j0(alpha*R0)
% for r discretized n times from 0 to R0
for i=1:n
    % go from 0 to R0
    r(i)=(i-1)*R0/(n-1);
    % calculate   1 - j0(alpha*r)/j0(alpha*R0)
    x(i)=unit(r(i),R0,freq);
    % unit gives real and imaginary part 
    % calculate the magnitude and phase
    m(i)=sqrt(real(x(i))^2 + imag(x(i))^2);
    if (~(real(x(i))==0)) 
	p(i)=atan(imag(x(i))/real(x(i))); 
    end;
end;

% show the mean of 1 - j0(a*r)/j0(a*R0) for changing r
a=mean(x);
mag=sqrt(real(a)^2 + imag(a)^2);
phase=atan(imag(a)/real(a));
out(1)=real(a);
out(2)=imag(a);
out(3)=mag;
out(4)=phase;
