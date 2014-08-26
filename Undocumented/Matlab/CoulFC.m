function F=CoulFC(x,R)

F=0.*x;

for i=1:numel(x)
F(i)=((R^3)/(x(i)^2))*sign(x(i));

if abs(x(i))<R
    F(i)=x(i);
end
end
