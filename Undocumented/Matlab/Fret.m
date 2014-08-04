function [Fret, tr, r, a, u, v] = Fret(pred1, pred2, t, acc)

r=zeros(numel(t),3);
a=zeros(numel(t),3);
u=zeros(numel(t),3);
v=zeros(numel(t),3);
Fret=zeros(numel(t),3);

rt = [(pred2.x - pred1.x), (pred2.y - pred1.y), (pred2.z - pred1.z)];

for i=2:numel(t)
    if t(i) - (sqrt(rt(i-1,:)*rt(i-1,:)')/PhysConst.c) > t(round(numel(t)/10))
        break
    end
end

beginsim = i;

tr = 0.*t;
itr = 0.*t;

v2 = [pred2.vx, pred2.vy, pred2.vz];

ax = diff(pred1.vx)./diff(t);
ay = diff(pred1.vy)./diff(t);
az = diff(pred1.vz)./diff(t);

disp('Start calculating tr')

options = optimoptions('fsolve','Display','off', 'TolX', acc);

for i=beginsim+1:numel(t)-1
    %time indication
    if i==round(0.05*(numel(t) - 2 - beginsim))
        disp('5% done')
    end
        if i==round(0.25*(numel(t) - 2 - beginsim))
        disp('25% done')
        end
        if i==round(0.50*(numel(t) - 2 - beginsim))
        disp('50% done')
        end
        if i==round(0.75*(numel(t) - 2 - beginsim))
        disp('75% done')
        end

    rtry = r(i-1,:);
    tr(i) = t(i) - ((sqrt(rtry*rtry'))/PhysConst.c);
    %multiply by a constant (scaling so a * r/c = 1 => a = c /r
    scale = PhysConst.c/((sqrt(rtry*rtry')));
    if scale > 1e20
        scale=1;
    end
    trsc = scale*tr(i);
    tsc = scale*t(i);
    rtrysc = scale*rtry;
    rtry2 = [pred2.x(i) - evalu(pred1.x, t, trsc/scale), pred2.y(i) - evalu(pred1.y, t, trsc/scale), pred2.z(i) - evalu(pred1.z, t, trsc/scale)];
    rtry2sc = scale*rtry2;
    %using scaling, tr -> c tr = 1 => c = 1/tr(i)
    n=0;
    while (((rtrysc - rtry2sc)*(rtrysc - rtry2sc)')/(rtrysc*rtrysc') > acc)
        n=n+1;
        rtrysc = rtry2sc;
        trsc = tsc - ((sqrt(rtrysc*rtrysc'))/PhysConst.c);
        rtry2sc = scale*[pred2.x(i) - evalu(pred1.x, t, trsc/scale), pred2.y(i) - evalu(pred1.y, t, trsc/scale), pred2.z(i) - evalu(pred1.z, t, trsc/scale)];
        if n>1e2
            rtry2sc=scale*fsolve(@(x)TrFind(x, pred1, pred2, t, i),rtrysc, options);
            break
            %try matlab solver
        end
    end
    tr(i) = trsc / scale;
    r(i,:) = rtry2sc./scale;
end

for i=beginsim+1:numel(t)-1

a(i,:) = [evalu(ax, t(i), tr(i)),evalu(ay, t(i), tr(i)),evalu(az, t(i), tr(i))];
v(i,:) = [evalu(pred1.vx, t(i), tr(i)),evalu(pred1.vy, t(i), tr(i)),evalu(pred1.vz, t(i), tr(i))];
u(i,:) = PhysConst.c.*r(i,:)./(sqrt(r(i,:)*r(i,:)')) - v(i,:);

Fret(i,:) =  (PhysConst.e^2 / (4 * pi * PhysConst.epsilon0)) * (sqrt(r(i,:)*r(i,:)')/((r(i,:)*u(i,:)').^3)) .* (((PhysConst.c^2 - v(i,:)*v(i,:)').*u(i,:) + ...
    cross(r(i,:), cross(u(i,:), a(i,:)))) + (1/PhysConst.c) .* ... 
    cross(v2(i,:), cross(r(i,:) ./sqrt(r(i,:)*r(i,:)'), (PhysConst.c^2 -  v(i,:)*v(i,:)').*u(i,:) + cross(r(i,:), cross(u(i,:), a(i,:))))));

    if tr(i) < 0
        tr(i)=0;
        Fret(i,:) = [0,0,0];
        a(i,:) = [0,0,0];
        v(i,:) = [0,0,0];
        u(i,:) = [0,0,0];
    end

if itr(i) == 1
    Fret(i,:) = [0, 0, 0];
end

end
