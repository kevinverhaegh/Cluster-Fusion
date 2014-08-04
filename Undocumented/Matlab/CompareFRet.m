function output = CompareFRet(trajec1, trajec2, acc, t)

rhelp=zeros(numel(t),3);
vhelp = zeros(numel(t),3);
FcInst = zeros(numel(t),3);

[FcRet21, tr21, r21, a21, u21, v21] = Fret(trajec1, trajec2, t, acc);
[FcRet12, tr12, r12, a12, u12, v12] = Fret(trajec2, trajec1, t, acc);

v2 = [trajec2.vx, trajec2.vy, trajec2.vz];

for i=1:numel(t)
    rhelp(i,:) = [trajec2.x(i) - trajec1.x(i), trajec2.y(i) - trajec1.y(i), trajec2.z(i) - trajec1.z(i)];
    vhelp(i,:) = [trajec1.vx(i), trajec1.vy(i),trajec1.vz(i)];
FcInst(i,:) =  (PhysConst.e^2 / (4 * pi * PhysConst.epsilon0))*(1/(rhelp(i,:)*rhelp(i,:)')) .* (rhelp(i,:)./sqrt(rhelp(i,:)*rhelp(i,:)') + 1/PhysConst.c^2 * ...
    cross(v2(i,:), cross(vhelp(i,:), rhelp(i,:)./sqrt(rhelp(i,:)*rhelp(i,:)'))));
end

FcInst12 = FcInst;
FcInst21 = - FcInst;

figure
plot(t,FcRet21(:,1),t,FcInst21(:,1), t, FcRet12(:,1),t,FcInst12(:,1))
xlabel('Time')
ylabel('Coulomb force x')
legend('Retarded 2->1', 'Instant 2->1', 'Retarded 1->2', 'Instant 1->2')
figure
plot(t,FcRet21(:,2),t,FcInst21(:,2), t, FcRet12(:,2),t,FcInst12(:,2))
xlabel('Time')
ylabel('Coulomb force y')
legend('Retarded 2->1', 'Instant 2->1', 'Retarded 1->2', 'Instant 1->2')
figure
plot(t,FcRet21(:,3),t,FcInst21(:,3), t, FcRet12(:,3),t,FcInst12(:,3))
xlabel('Time')
ylabel('Coulomb force z')
legend('Retarded 2->1', 'Instant 2->1', 'Retarded 1->2', 'Instant 1->2')

output=struct('Fret21', FcRet21, 'Fdir12', FcInst, 'tr21', tr21, 't', t, 'Fret12', FcRet12, 'tr12', tr12);
