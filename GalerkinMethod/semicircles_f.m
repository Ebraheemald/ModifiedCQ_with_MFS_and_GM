function p = semicircles_f(M)

delta = (2*pi)/(M);  % the circumference is 2pi for the unit circle
a1 = (3*pi)/2;   b1 = (pi*5)/2;
a2 = pi/2;       b2 = (3*pi)/2;

% == upper Right semicirlce ==
s1 = a1:delta:b1;
p1 = exp(1i*s1);      

% == Upper semicirlce ========
d = .25;
s2 = a2:4*delta:b2;
p2a = exp(1i*s2);
p2 = p2a*d; p2 = p2+1i-d*1i; 

% == Left semicirlce =========
s3 = a1:2*delta:b1;
p3a = 0.5*exp(1i*s3);
p3 = flip(p3a);

% == Lower semicirlce ========
p4 = p2a*d; p4 = p4-1i+d*1i;   

% == Total domain ============
p = [p1(1:end-1) p2(1:end-1) p3(1:end-1) p4(1:end-1)];
p = -0.5+.125+p;