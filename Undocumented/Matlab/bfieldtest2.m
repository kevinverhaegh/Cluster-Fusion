%Calculates the magnetic field using finite line segments. Assumption is
%that z=0. 

%make x, y-vector, discretizing positions
x=-12:1.1:12;
y=-12:1.1:12;

%general equation
%A = (mu0*I)/(4*pi) * Log((z2 + Sqrt(z2^2 + s^2))/(z1 + Sqrt(z1^2 + s^2))) * Direction

%From each coordinate transformation follows matrices for z2, s, Direction

%constants
mu0 = 4*pi*1e-7;
I=1e8;

%the direction of the vector potential of each line segment
Direction = [0 1; 0 -1; 1 0; -1 0; 0.5*sqrt(2) 0.5*sqrt(2); -0.5*sqrt(2) -0.5*sqrt(2); 0.5*sqrt(2) -0.5*sqrt(2); -0.5*sqrt(2) 0.5*sqrt(2)];

%Setting matrices zero
Ax = zeros(numel(x), numel(y));
Ay = zeros(numel(x), numel(y));
z2=zeros(numel(Direction),1);
z1=zeros(numel(Direction),1);
s=zeros(numel(Direction),1);

%looping over all positions
for i=1:1:numel(x)
for j=1:1:numel(y)

%the position of the general z2
z2 = [y(j)+5, y(j)+5, x(i)+5, x(i)+5, 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)-7.5) + 0.5*sqrt(2)*5, 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)+7.5) + 0.5*sqrt(2)*5, - 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)-7.5) + 0.5*sqrt(2)*5, - 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)+7.5) + 0.5*sqrt(2)*5];

%the position of the general z1
z1 = [y(j)-5, y(j)-5, x(i)-5, x(i)-5, 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)-7.5) - 0.5*sqrt(2)*5, 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)+7.5) - 0.5*sqrt(2)*5, - 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)-7.5) - 0.5*sqrt(2)*5, - 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)+7.5) - 0.5*sqrt(2)*5];

%the position of the general s
s = [x(i)+10, x(i)-10, y(j)-10, y(j)+10, 0.5*sqrt(2)*(x(i)+7.5) - 0.5*sqrt(2)*(y(j)-7.5), 0.5*sqrt(2)*(x(i)-7.5) - 0.5*sqrt(2)*(y(j)+7.5), 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)-7.5), 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)+7.5)];
    
%looping over each line segment
for k=1:1:8

%adding all vector potentials of each line segment
Ax(i,j) = Ax(i,j) + (mu0*I)/(4*pi) * log((z2(k) + sqrt(z2(k)^2 + s(k)^2))/(z1(k) + sqrt(z1(k)^2 + s(k)^2))) * Direction(k,1);
Ay(i,j) = Ay(i,j) + (mu0*I)/(4*pi) * log((z2(k) + sqrt(z2(k)^2 + s(k)^2))/(z1(k) + sqrt(z1(k)^2 + s(k)^2))) * Direction(k,2);
end
end
end

%draws vector potential
quiver(x,y, Ax , Ay )

%calculates curl of the vector potential: B = curl(A) = nabla x A
[curlz,cav]= curl(y,x,Ay,Ax);
%curlz is B in the z direction
figure
%plots magnetic flux density
surf(x,y,curlz, 'EdgeColor','none','LineStyle','none')