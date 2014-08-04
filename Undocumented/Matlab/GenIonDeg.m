function [Deg, ElC, Zt, rm,trapEtab]=GenIonDeg(output)
Zt=zeros(numel(output.t),1);
rm=zeros(numel(output.t),1);
ElC=zeros(numel(output.t),1);
trapEtab=zeros(numel(output.t),numel(output.x{1,1}(1,:)),3);
for i=1:numel(output.t)
    for j=1:numel(output.GroupsClass(1,:))
        if output.GroupsClass(i,j)==3
            Zt(i)=Zt(i)+output.AMUZcombo(output.GroupsClass(i,j),2);
                if abs(output.x{output.GroupsClass(i,j),1}(i,j)) > rm(i)
                rm(i) = abs(output.x{output.GroupsClass(i,j),1}(i,j));
                end    
        end
    end
    for j=1:numel(output.GroupsClass(2,:))
    if output.GroupsClass(i,j)==1 && sqrt(output.x{1,1}(i,j)^2 + output.y{1,1}(i,j)^2 + output.z{1,1}(i,j)^2)<rm(i)
        %trapped electron detected
        trapEtab(i,j,1) = output.x{1,1}(i,j);
        trapEtab(i,j,2) = output.y{1,1}(i,j);
        trapEtab(i,j,3) = output.z{1,1}(i,j);
        ElC(i)=ElC(i)+1;
    end
    end
    disp(100 * i/numel(output.t))
end
plot(output.t, 1-ElC./Zt)
Deg = 1-ElC./Zt;
