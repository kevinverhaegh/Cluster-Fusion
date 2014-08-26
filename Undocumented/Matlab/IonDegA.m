function Deg=IonDegA(output,j)

IonDegH=0;
EkH=0;
Rm=sqrt(max(output.x{3,1}(j,:).^2 + output.y{3,1}(j,:).^2 + output.z{3,1}(j,:).^2));
for i=1:numel(output.x{1,1}(1,:))
    if output.x{1,1}(end,i)~= 0
        if sqrt(output.x{1,1}(j,i)^2 + output.y{1,1}(j,i)^2 + output.z{1,1}(j,i)^2) < Rm
            IonDegH=IonDegH+1;
            EkH=EkH+PhysConst.me*PhysConst.c^2 * ((1 / sqrt( 1 - (output.vx{1,1}(j,i)^2 + output.vy{1,1}(j,i)^2 + output.vz{1,1}(j,i)^2)))-1);
        end
    end
end
Deg=1- IonDegH/output.NumberInGroup(j,1);