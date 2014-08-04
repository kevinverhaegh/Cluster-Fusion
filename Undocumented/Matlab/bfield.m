%make x, y-vector
x=-12:1:12;
y=-12:1:12;

%general equation
%A = (mu0*I)/(4*pi) * Log((z2 + Sqrt(z2^2 + s^2))/(z1 + Sqrt(z1^2 + s^2))) * Direction

%From each coordinate transformation follows matrices for z2, s, Direction

%constants
mu0 = 4*pi*1e-7;
I=1e8;

Direction = [0 1; 0 -1; 1 0; -1 0; 0.5*sqrt(2) 0.5*sqrt(2); -0.5*sqrt(2) -0.5*sqrt(2); 0.5*sqrt(2) -0.5*sqrt(2); -0.5*sqrt(2) 0.5*sqrt(2)];

A = zeros(numel(x), numel(y), 2);

z2=zeros(8,1);
z1=zeros(8,1);
s=zeros(8,1);

for i=1:1:numel(x)
for j=1:1:numel(y)

z2 = [y(j)-5, y(j)+5, x(i)-5, x(i)-5, 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)-7.5) - 0.5*sqrt(2)*5, 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)+7.5) - 0.5*sqrt(2)*5, - 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)-7.5) - 0.5*sqrt(2)*5, - 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)+7.5) - 0.5*sqrt(2)*5];
    
z1 = [y(j)+5, y(j)-5, x(i)+5, x(i)+5, 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)-7.5) + 0.5*sqrt(2)*5, 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)+7.5) + 0.5*sqrt(2)*5, - 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)-7.5) + 0.5*sqrt(2)*5, - 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)+7.5) + 0.5*sqrt(2)*5];
    
s = [x(i)+10, x(i)-10, y(j)+10, y(j)-10, 0.5*sqrt(2)*(x(i)+7.5) - 0.5*sqrt(2)*(y(j)-7.5), 0.5*sqrt(2)*(x(i)-7.5) - 0.5*sqrt(2)*(y(j)+7.5), 0.5*sqrt(2)*(x(i)-7.5) + 0.5*sqrt(2)*(y(j)-7.5), 0.5*sqrt(2)*(x(i)+7.5) + 0.5*sqrt(2)*(y(j)+7.5)];
    
for k=1:1:numel(s)

A(i,j,:) = A(i,j,:) + (mu0*I)/(4*pi) * log((z2(k) + sqrt(z2(k)^2 + s(k)^2))/(z1(k) + sqrt(z1(k)^2 + s(k)^2))) * Direction(k);

end
end
end

%Remove NaNs and Infs @ the boundaries

for i=1:1:numel(x)
for j=1:1:numel(y)
for k=1:1:2
    
    if isnan(A(i,j,k)) == 1 
        A(i,j,k) = 0;
    end
    if isinf(A(i,j,k)) ==1
        A(i,j,k) = 0;
    end
end
end
end

quiver(x,y, A(:,:,1) , A(:,:,2) )