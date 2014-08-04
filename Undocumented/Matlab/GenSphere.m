function [x,y,z]=GenSphere(N,R)

x=zeros(N,1);
y=0.*x;
z=0.*x;

help =randsphere(N-1,3,R);
x(2:N)=help(:,1);
y(2:N)=help(:,2);
z(2:N)=help(:,3);

    