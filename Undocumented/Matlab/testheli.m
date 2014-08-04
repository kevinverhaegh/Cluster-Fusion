function output =testheli(env, delta, tb, te)

t = tb:((te-tb)/(1e5)):te;

pred1 = ElectronField('GeneralTravel', 0,0,0,env,0,0,0,t);
dz = delta*(2*pi*PhysConst.c/env.omega);
pred2 = ElectronField('GeneralTravel', 0,0,dz,env,0,0,0,t);

frac = 0.7;

R = 1e-6;
omega = frac*PhysConst.c/R;

fracH = 0.7;
vH = fracH.*PhysConst.c;

pred1.x = 0.*t;
pred1.y = 0.*t;
pred1.z = vH.*t;
pred1.vz = vH.*t./t;
pred1.vy = 0.*t;
pred1.vx = 0.*t;

pred2.z = vH.*t;
pred2.y = R.*sin(omega.*t);
pred2.x = R.*cos(omega.*t);
pred2.vz = vH.*t./t;
pred2.vy = R.*omega.*cos(omega.*t);
pred2.vx = - R.*omega.*sin(omega.*t);

switchparam = 1;

if switchparam == 1
    dum1 = pred1;
    dum2 = pred2;
    pred1 = dum2;
    pred2 = dum1;
end

rt = [(pred2.x - pred1.x)', (pred2.y - pred1.y)', (pred2.z - pred1.z)'];

for i=2:numel(t)
    if t(i) - (sqrt(rt(i-1,:)*rt(i-1,:)')/PhysConst.c) > t(1)
        break
    end
end

beginsim = i;

%adjust t to t adaptive
t0 = tb:((te-tb)/(1e5)):t(beginsim)-((te-tb)/(1e5));
t1 = t(beginsim):1e-20:t(beginsim)+1e-16 -1e-20;
t2 = t(beginsim)+1e-16:((te-tb)/(1e5)):te;
clear t
t = horzcat(t0, t1, t2);

clear pred1
pred1 = ElectronField('GeneralTravel', 0,0,0,env,0,0,0,t);
dz = delta*(2*pi*PhysConst.c/env.omega);
clear pred2
pred2 = ElectronField('GeneralTravel', 0,0,dz,env,0,0,0,t);

r=zeros(numel(t),3);
a=zeros(numel(t),3);
u=zeros(numel(t),3);
v=zeros(numel(t),3);
FcRet=zeros(numel(t),3);
FcInst = zeros(numel(t),3);
rhelp = zeros(numel(t),3);
vhelp=zeros(numel(t),3);

%vt = 0.9*PhysConst.c;

pred1.x = 0.*t;
pred1.y = 0.*t;
pred1.z = vH.*t;
pred1.vx = 0.*t;
pred1.vy = 0.*t;
pred1.vz = vH.*t./t;

pred2.z = vH.*t;
pred2.y = R.*sin(omega.*t);
pred2.x = R.*cos(omega.*t);
pred2.vz = vH.*t./t;
pred2.vy = R.*omega.*cos(omega.*t);
pred2.vx = - R.*omega.*sin(omega.*t);

if switchparam == 1
    dum1 = pred1;
    dum2 = pred2;
    pred1 = dum2;
    pred2 = dum1;
end

tr = 0.*t;
itr = 0.*t;

v2 = [pred2.vx', pred2.vy', pred2.vz'];

ax = diff(pred1.vx)./diff(t);
ay = diff(pred1.vy)./diff(t);
az = diff(pred1.vz)./diff(t);

for i=beginsim+1:numel(t)
    rtry = r(i-1,:);
    tr(i) = t(i) - ((sqrt(rtry*rtry'))/PhysConst.c);
    [~,itr(i)] = min(abs(tr(i) - t));
    rtry2 = [pred2.x(i) - pred1.x(itr(i)), pred2.y(i) - pred1.y(itr(i)), pred2.z(i) - pred1.z(itr(i))];
    n=0;
    while (((rtry - rtry2)*(rtry - rtry2)')/(rtry*rtry') > 1e-10)
        n=n+1;
        rtry = rtry2;
        tr(i) = t(i) - ((sqrt(rtry*rtry'))/PhysConst.c);
        [~,itr(i)] = min(abs(tr(i) - t));
        rtry2 = [pred2.x(i) - pred1.x(itr(i)), pred2.y(i) - pred1.y(itr(i)), pred2.z(i) - pred1.z(itr(i))];
        if n>1e5
            disp('Cannot converge')
            break
        end
    end
    if n>1e5
        break
    end
    
    r(i,:) = rtry2;

a(i,:) = [ax(itr(i)),ay(itr(i)),az(itr(i))];
v(i,:) = [pred1.vx(itr(i)), pred1.vy(itr(i)), pred1.vz(itr(i))];
u(i,:) = PhysConst.c.*r(i,:)./(sqrt(r(i,:)*r(i,:)')) - v(i,:);

FcRet(i,:) =  (PhysConst.e^2 / (4 * pi * PhysConst.epsilon0)) * (sqrt(r(i,:)*r(i,:)')/((r(i,:)*u(i,:)').^3)) .* (((PhysConst.c^2 - v(i,:)*v(i,:)').*u(i,:) + ...
    cross(r(i,:), cross(u(i,:), a(i-1,:)))) + (1/PhysConst.c) .* ... 
    cross(v2(i,:), cross(r(i,:) ./sqrt(r(i,:)*r(i,:)'), (PhysConst.c^2 -  v(i,:)*v(i,:)').*u(i,:) + cross(r(i,:), cross(u(i,:), a(i,:))))));

    if tr(i) < 0
        tr(i)=0;
        FcRet(i,:) = [0,0,0];
        a(i,:) = [0,0,0];
        v(i,:) = [0,0,0];
        u(i,:) = [0,0,0];
    end

if itr(i) == 1
    FcRet(i,:) = [0, 0, 0];
end

%disp(i/numel(t))
end

for i=1:numel(t)
    rhelp(i,:) = [pred2.x(i) - pred1.x(i), pred2.y(i) - pred1.y(i), pred2.z(i) - pred1.z(i)];
    vhelp(i,:) = [pred1.vx(i), pred1.vy(i),pred1.vz(i)];
FcInst(i,:) =  (PhysConst.e^2 / (4 * pi * PhysConst.epsilon0))*(1/(rhelp(i,:)*rhelp(i,:)')) .* (rhelp(i,:)./sqrt(rhelp(i,:)*rhelp(i,:)') + 1/PhysConst.c^2 * ...
    cross(v2(i,:), cross(vhelp(i,:), rhelp(i,:)./sqrt(rhelp(i,:)*rhelp(i,:)'))));
end

figure
plot(t,FcRet(:,1),t,FcInst(:,1))
xlabel('Time')
ylabel('Coulomb force x')
legend('Retarded potentials', 'Instant Coulomb')
figure
plot(t,FcRet(:,2),t,FcInst(:,2))
xlabel('Time')
ylabel('Coulomb force y')
legend('Retarded potentials', 'Instant Coulomb')
figure
plot(t,FcRet(:,3),t,FcInst(:,3))
xlabel('Time')
ylabel('Coulomb force z')
legend('Retarded potentials', 'Instant Coulomb')

output=struct('Fret', FcRet, 'Fdir', FcInst, 'tr', tr, 't', t, 'x',pred1.x, 'v', pred1.vx);
