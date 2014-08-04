function dy = differ(t,y)
M=5;
m=1;
J=1e-2;
l=0.1;
g=9.81;
xi = J + m*(l^2) - (m^2 * l^2)/(M+m);
a = (m*g*l)/xi
b = ((m*l)/(m+M))/xi
Kd=0.5*b;
Kp=1.1*a;

dy = zeros(2,1);
%u=(1e-5)*sin(2*pi*t);
u=(Kp/b)*y(1)+(Kd/b)*y(2);
dy(1) = y(2);
dy(2) = (m*g*l*sin(y(1)) - ((m*l)/M)*u*cos(y(1)))/(J+m*(l^2));
%dy(2) = (-((m*l)/(M+m))*u + m*l*g*y(1))/(m*(l^2) + J - ((m*l)^2)/(M+m));


end

