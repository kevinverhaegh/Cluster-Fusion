function PlotField(Fx,Fy,Fz,x,y,z,tind)
close all
figure
quiver3(x',(zeros(numel(x),1)+1).*y(1),(zeros(numel(x),1)+1).*z(1),Fx(:,1,1,tind),Fy(:,1,1,tind),Fz(:,1,1,tind))
hold
for i=2:numel(y)
    for j=2:numel(z)
        quiver3(x',(zeros(numel(x),1)+1).*y(i),(zeros(numel(x),1)+1).*z(j),Fx(:,i,j,tind),Fy(:,i,j,tind),Fz(:,i,j,tind))
    end
end
    