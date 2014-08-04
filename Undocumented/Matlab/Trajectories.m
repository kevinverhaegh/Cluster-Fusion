function output = Trajectories(info)
size = numel(info);
N = info{1,1}.param.nps;
groups = numel(info{1,1}.data.x)/N;
groupsID = zeros(groups,1);
for i=1:groups
    groupsID(i) = info{1,1}.data.q(1+N*(i-1))/PhysConst.e;
end

x = zeros(size, N, groups);
y = zeros(size, N, groups);
z = zeros(size, N, groups);
vx = zeros(size, N, groups);
vy = zeros(size, N, groups);
vz = zeros(size, N, groups);
fEx = zeros(size,N,groups);
fEy = zeros(size,N,groups);
fEz = zeros(size,N,groups);
fBx = zeros(size,N,groups);
fBy = zeros(size,N,groups);
fBz = zeros(size,N,groups);
t = zeros(size,1);

for i=1:size
    t(i) = info{i,1}.param.time;
        for j=1:groups
            x(i,:,j) = info{i,1}.data.x((j-1)*N+1:j*N);
            y(i,:,j) = info{i,1}.data.y((j-1)*N+1:j*N);
            z(i,:,j) = info{i,1}.data.z((j-1)*N+1:j*N);
            vx(i,:,j) = info{i,1}.data.Bx((j-1)*N+1:j*N).*PhysConst.c;
            vy(i,:,j) = info{i,1}.data.By((j-1)*N+1:j*N).*PhysConst.c;
            vz(i,:,j) = info{i,1}.data.Bz((j-1)*N+1:j*N).*PhysConst.c;
            fEx(i,:,j) = info{i,1}.data.fEx((j-1)*N+1:j*N);
            fEy(i,:,j) = info{i,1}.data.fEy((j-1)*N+1:j*N);
            fEz(i,:,j) = info{i,1}.data.fEz((j-1)*N+1:j*N);
            fBx(i,:,j) = info{i,1}.data.fBx((j-1)*N+1:j*N);
            fBy(i,:,j) = info{i,1}.data.fBy((j-1)*N+1:j*N);
            fBz(i,:,j) = info{i,1}.data.fBz((j-1)*N+1:j*N);
        end
end

output.x=x;
output.y=y;
output.z=z;
output.vx=vx;
output.vy=vy;
output.vz=vz;
output.fEx=fEx;
output.fEy=fEy;
output.fEz=fEz;
output.fBx=fBx;
output.fBy=fBy;
output.fBz=fBz;
output.t = t;

         
    
    