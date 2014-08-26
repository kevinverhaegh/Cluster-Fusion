E=1e9:1e9:5e10;
x=zeros(numel(E),100);
for i=1:numel(E)
    for j=1:100
    x0=rand(1)*1e-9 - 5e-10;
    xx=Oscillatory(E(i),x0);
    x(i,j)=xx(end,1);
    end
end

