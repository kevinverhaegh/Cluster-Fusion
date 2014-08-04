function output=Etot(input)
vx=input.vx;
vy=input.vy;
vz=input.vz;
x=input.x;
y=input.y;
z=input.z;
GroupID=input.GroupID;
q=0.*vx;
for i=1:numel(GroupsID)
q(:,:,i)=(GroupsID(i) + zeros(numel(vx(:,1,1)),numel(vx(1,:,1)))).*PhysConst.e;
end
mdum=input.m;
m=0.*vx;
for i=1:numel(GroupsID)
m(:,:,i)=(m(i) + zeros(numel(vx(:,1,1)),numel(vx(1,:,1))));
end
gamma = 1/sqrt(1 - (vx.^2 + vy.^2 + vz.^2)./PhysConst.c^2);
Ekin = (gamma - 1).*m.*(PhysConst.c^2);
for i=1:numel(groupsID)