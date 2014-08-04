function [Pred, Info,Analysis] = SpherePred(type,env,t,N,R)

[x,y,z]=GenSphere(N,R);

Pred=cell(N,1);

xp=zeros(numel(t),N);
yp=zeros(numel(t),N);
zp=zeros(numel(t),N);

vx=zeros(numel(t),N);
vy=zeros(numel(t),N);
vz=zeros(numel(t),N);

for i=1:N
    
    Pred{i,1}=ElectronField(type, x(i),y(i),z(i),env,0,0,0,t);
    xp(:,i)=Pred{i,1}.x;
    yp(:,i)=Pred{i,1}.y;
    zp(:,i)=Pred{i,1}.z;
    vx(:,i) = Pred{i,1}.vx;
    vy(:,i) = Pred{i,1}.vy;
    vz(:,i)=Pred{i,1}.vz;
end

xg = 0.*t;
yg = 0.*t;
zg = 0.*t;
sxg = 0.*t;
syg = 0.*t;
szg = 0.*t;
vxg = 0.*t;
vyg = 0.*t;
vzg = 0.*t;
svxg = 0.*t;
svyg = 0.*t;
svzg = 0.*t;

for i=1:numel(t)
    xg(i) = mean(xp(i,:));
    sxg(i) = std(xp(i,:));
    yg(i) = mean(yp(i,:));
    syg(i) = std(yp(i,:));
    zg(i) = mean(zp(i,:));
    szg(i) = std(zp(i,:));
    vxg(i) = mean(vx(i,:));
    svxg(i) = std(vx(i,:));
    vyg(i) = mean(vy(i,:));
    svyg(i) = std(vy(i,:));
    vzg(i) = mean(vz(i,:));
    svzg(i) = std(vz(i,:));
end

Info=struct('xp',xp,'yp',yp, 'zp',zp,'vx',vx,'vy',vy,'vz',vz);

Analysis=struct('xg',xg,'sxg',sxg,'yg',yg,'syg',syg,'zg',zg,'szg',szg, 'vxg',vxg,'svxg',svxg,'vyg',vyg,'svyg',svyg,'vzg',vzg,'svzg',svzg);


