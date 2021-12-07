function V = slp_mod(ss1,ss2,r,dt,t)

mij=floor(r/dt); tm = t(mij+1);
V= (ss1.^mij).*exp(ss2*(tm-r))*(1/(2*pi)).*besselk(0,ss2*r,1);