function [ColMat, VrelMat, ColCount]=ColInv1(input,tma, tmi, coldefR, coldefV)

[~,indmi]=min(abs(tmi-input.t));
[~,indm]=min(abs(tma-input.t));
ColMat = cell(indm-indmi,1);
VrelMat = cell(indm-indmi,1);
for i=indmi:indm
    ColMat{i,1}=sparse(zeros(max(max(input.NumberInGroup)),max(max(input.NumberInGroup))));
    VrelMat{i,1}=sparse(zeros(max(max(input.NumberInGroup)),max(max(input.NumberInGroup))));
end
ColCount = zeros(indm-indmi,1);

N=max(max(input.NumberInGroup));

try
for i=indmi:indm
    xjS=full(input.x{1,1}(i,1:N));
    yjS=full(input.y{1,1}(i,1:N));
    zjS=full(input.z{1,1}(i,1:N));
    vxjS=full(input.vx{1,1}(i,1:N));
    vyjS=full(input.vy{1,1}(i,1:N));
    vzjS=full(input.vz{1,1}(i,1:N));  
    for j=1:N
        r=sqrt((xjS(j)-xjS).^2 + (yjS(j)-yjS).^2 + (zjS(j)-zjS).^2);
        vrel=sqrt((vxjS(j)-vxjS).^2 + (vyjS(j)-vyjS).^2 + (vzjS(j)-vzjS).^2);
        
        Bool1 = r<coldefR;
        Bool2 = vrel>coldefV;
        
        Bool=Bool1.*Bool2;
        
        ColMat{i,1}(j,:)=Bool;
        VrelMat{i,1}(j,:)=Bool.*vrel;
    end
        disp(i/(indm-indmi));
    ColCount(i)=sum(sum(ColMat{i,1}));
end

catch
    disp('Fail');
end

        