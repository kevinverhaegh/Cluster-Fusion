syms p(t)

M=5;
m=1;
l=0.1;
J=1e-2;
g=9.81;

u = (1e-5)*sin(2*pi*t);
phin = 0;

sdiff = u/(M+m) - ((m*l)/(M+m))*diff(p,2)*cos(p) + ((m*l)/(M+m))*(diff(p)^2)*sin(p);

S=dsolve(J*diff(p,2) == - m*sdiff*l*cos(p) - m*(l^2)*diff(p,2) + m*g*sin(p)*l, p(0)==phin);


