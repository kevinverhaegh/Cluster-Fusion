
%symbol toolkit
syms q1 q2 epsilon0 c x1(t) y1(t) z1(t) x2(t) y2(t) z2(t)  

%Constant
C = q1*q2 / (4*pi*epsilon0);

%Position particle 1
xv1 = [x1, y1, z1];

%Position particle 2
xv2 = [x2, y2, z2];

%difference vectors, u vectors
rv12 = xv2 - xv1;
rv21 = - rv12;
xv1d = diff(xv1, t);
xv2d = diff(xv2, t);
xv1dd = diff(xv1d, t);
xv2dd = diff(xv2d, t);
u12 = (c .* rv12 ./ (sqrt(rv12*rv12'))) - xv1d;
u21 = (c .* rv21 ./ (sqrt(rv21*rv21'))) - xv2d;

%Retarded Coulomb forces (12 means from particle 1 working on 2, etc.)
Fc12 = C * (sqrt(rv12(tr)*rv12(tr)')/((rv12(tr).*u)^3)) * ...
    (((c^2 - xv1d*xv1d') .* u + cross(rv12, cross(u, xv1dd))) ...
    + (1/c)*cross(xv2d, ...
    cross(rv12, ...
    (c^2 - xv1d*xv1d').*u + cross(rv12, cross(u, xv1dd)))));
Fc21 = C * (sqrt(rv12(tr)*rv12(tr)')/((rv21(tr).*u)^3)) * ...
    (((c^2 - xv2d*xv2d') .* u + cross(rv21, cross(u, xv2dd))) ...
    + (1/c)*cross(xv1d, ...
    cross(rv21, ...
    (c^2 - xv2d*xv2d').*u + cross(rv21, cross(u, xv2dd)))));

%External field forces
Ff1 = q1*E(xv1,t) + q1 * cross(xv1d,B(xv1,t));
Ff2 = q2*E(xv2,t) + q2 * cross(xv2d,B(xv2,t));

%Lorentz factor
gamma1 = 1/sqrt(1 - ((xv1d*xv1d')/c));
gamma2 = 1/sqrt(1 - ((xv2d*xv2d')/c));

%momentum derivative
dmom1 = gamma1*m1.*(xv1d*xv1dd' * gamma1^2 .* xv1d + xv1dd);
dmom2 = gamma2*m2.*(xv2d*xv2dd' * gamma2^2 .* xv2d + xv2dd);