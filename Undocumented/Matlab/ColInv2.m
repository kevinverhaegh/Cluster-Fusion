function [ColMat, VrelMat]=ColInv2(input,tma, tmi, coldefR, coldefV)

[~,indmi]=min(abs(tmi-input.t));
[~,indm]=min(abs(tma-input.t));
ColMat = zeros(indmi-indm, max(max(input.NumberInGroup)),max(max(input.NumberInGroup)));
VrelMat = zeros(indmi-indm, max(max(input.NumberInGroup)),max(max(input.NumberInGroup)));
ColCount = zeros(idmi-indm,1);

N=max(max(input.NumberInGroup));

indNZv=0;

try
for i=indmi:indm
    xjS=full(input.x{1,1}(i,N+1:end));
    yjS=full(input.y{1,1}(i,N+1:end));
    zjS=full(input.z{1,1}(i,N+1:end));
    vxjS=full(input.vx{1,1}(i,N+1:end));
    vyjS=full(input.vy{1,1}(i,N+1:end));
    vzjS=full(input.vz{1,1}(i,N+1:end));
    
    xjS2=full(input.x{2,1}(i,1:N));
    yjS2=full(input.y{2,1}(i,1:N));
    zjS2=full(input.z{2,1}(i,1:N));
    vxjS2=full(input.vx{2,1}(i,1:N));
    vyjS2=full(input.vy{2,1}(i,1:N));
    vzjS2=full(input.vz{2,1}(i,1:N));
    
    for j=1:N
        xj=xjS(j);
        yj=yjS(j);
        zj=zjS(j);
        vxj=vxjS(j);
        vyj=vyjS(j);
        vzj=vzjS(j);
        
        r=sqrt((xj-xjS2).^2 + (yj-yjS2).^2 + (zj-zjS2).^2);
        vrel=sqrt((vxj-vxjS2).^2 + (vyj-vyjS2).^2 + (vzj-vzjS2).^2);
        
        Bool1 = r<coldefR;
        Bool2 = vrel>coldefV;
        
        Bool=Bool1.*Bool2;
        
        ColMat(i,j,:)=Bool;
        VrelMat(i,j,:)=Bool.*vrel;
    end
    disp(i/(indm-indmi));
    ColCount(i)=sum(sum(ColMat(i,:,:)));
end

catch
    disp('Fail');
end

        