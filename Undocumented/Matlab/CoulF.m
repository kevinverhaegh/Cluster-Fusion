function F=CoulF(x,CutOff)

F=0.*x;

for i=1:numel(x)
F(i)=sign(x(i))*(CutOff/abs(CutOff))*(1/CutOff^2);

if abs(x(i))>CutOff
    F(i)=(sign(x(i)))*(1/x(i)^2);
end
end

