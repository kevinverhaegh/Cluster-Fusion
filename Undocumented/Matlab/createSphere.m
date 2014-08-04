function [actSize,x,y,z]=createSphere(N)
SphereToCube = 4/3 * pi * (0.5^3);
Nupgr=ceil(N/SphereToCube);

Nlupgr = ceil(Nupgr^(1/3));
xC = -1:2 * 1/Nlupgr:1;
yC = xC;
zC = xC;

counterr=0;
CubeInitial = zeros(numel(xC)^3, 3);
for i=1:numel(xC)
    for j=1:numel(xC)
        for k=1:numel(xC)
            counterr=counterr+1;
            CubeInitial(counterr,1)=xC(i);
            CubeInitial(counterr,2)=yC(j);
            CubeInitial(counterr,3)=zC(k);
        end
    end
end

SphereInitial=0.*CubeInitial;
counter = 0;

for i=1:numel(CubeInitial(:,1))
    if sqrt(CubeInitial(i,1).^2 + CubeInitial(i,2).^2 + CubeInitial(i,3).^2) < 1
        counter = counter + 1;
        SphereInitial(counter,1)=CubeInitial(i,1);
        SphereInitial(counter,2)=CubeInitial(i,2);
        SphereInitial(counter,3)=CubeInitial(i,3);
    end
end

actSize=counter;
x=SphereInitial(1:counter,1);
y=SphereInitial(1:counter,2);
z=SphereInitial(1:counter,3);
