function output=ShockShellModelNAG(rho, R0, Z, M, t, r0t)
global taui Qi s0i
N=trapz(r0t, rho.*.4.*pi.*r0t.^2);

disp('Amount of ions:')
disp(N)

densavg=N/(4/3 * pi * R0^3) ; 

t0 = (1/(N*PhysConst.e*Z))*sqrt(M*(R0^3)/N);
%t0 = sqrt(M*PhysConst.epsilon0/(densavg * Z^2 * PhysConst.e^2));

disp('CE time scale:')
disp(t0)

s0=r0t./R0;
tau = t./t0;

Q=cumtrapz(r0t, 4*pi*rho.*r0t.^2)/N;

rvecN=zeros(numel(rho), numel(t));
rhoN=zeros(numel(rho),numel(t));

try
fun=@(x,Qi,s0i,taui) (sqrt(x*(x-1)) + log(sqrt(x) + sqrt(x-1)) - taui*sqrt(2*Qi/(3*s0i^3)));
for i=1:numel(tau)
    taui=tau(i);
    for j=1:numel(s0)
        Qi=Q(j);
        s0i=s0(j);
        rvecN(j,i)=r0t(j)*real(c05nb(1,'f_solve_obj_NAG', 1+1e-6));
        %specfun = @(x) fun(x, Qi, s0i, taui) + (x>1);
        %opts = optimset('Diagnostics','off', 'Display','off', 'TolFun', 1e-9, 'TolX', 1e-9);
        %rvecN(j,i)=r0t(j)*real(fsolve(specfun, 1+1e-6,opts));
        rhoN(j,i) = real((rvecN(j,i)^-2)*rho(j));
    end
end


plot(rvecN')

vvecN=diff(rvecN')/(t(2)-t(1));
%each virtual particle attributes to a shell, which represents a number of
%particles, calculated here

NumPVir=cumtrapz([0 r0t], 4*pi*[rho(1) rho].*[0 r0t].^2);
NumPVir=[NumPVir(1) diff(NumPVir)];

cross=0;
totcross=0;
%search for collisions, a collision is an intersection of the various lines
%of various radia. 
for i=1:numel(s0)-1
    for j=i+1:numel(s0)
        icross = intersection(rvecN(i,:), rvecN(j,:));
        if icross~=0
        totcross=totcross+1;
        end
    end
end

icrossvec=zeros(totcross,1);
tcross=0.*icrossvec;
v1cross=0.*icrossvec;
v2cross=0.*icrossvec;
rcross=zeros(totcross,2);

for i=1:numel(s0)-1
    for j=i+1:numel(s0)
        icross = intersection(rvecN(i,:), rvecN(j,:));
        if icross~=0
            cross=cross+1;
            icrossvec(cross)=round(icross);
            rcross(cross, :)=[i, j];
            tcross(cross)=findv(t,icross);
            v1cross(cross)=findv(vvecN(:,i), icross);
            v2cross(cross)=findv(vvecN(:,j),icross);
        end
    end
end

%calculate fusion reactions
%load DD crossections
load('DD_crossection.mat');
reac=0.*v1cross;
Erel = 0.5 * M * abs(v1cross - v2cross).^2 ./ PhysConst.e;
for i=1:numel(reac)
if Erel(i)<20
reac(i)=0;
else
[~,r]=min(abs(Erel(i)-EnDDsp));
reac(i)=DDsp(r)*(10^-28)*abs(v1cross(i)-v2cross(i));
end
end

V1=reac.*0;
V2=reac.*0;
n1=0.*V1;
n2=0.*V2;
tdur=0.*V1;
FusY=0.*V1;
for i=1:numel(reac)
    if rcross(i,1)~=numel(rvecN(:,1)) && rcross(i,2)~=numel(rvecN(:,1))
    V1(i)=4*pi*(rvecN(rcross(i,1), icrossvec(i)).^2) * (rvecN(rcross(i,1)+1, icrossvec(i)) - rvecN(rcross(i,1), icrossvec(i)));
    else if rcross(i,1)~=1
        V1(i)=4*pi*(rvecN(rcross(i,1), icrossvec(i)).^2) * (rvecN(rcross(i,1), icrossvec(i)) - rvecN(rcross(i,1)-1, icrossvec(i)));
        end
    end
    
    if rcross(i,1)~=numel(rvecN(:,2)) && rcross(i,2)~=numel(rvecN(:,1))
    V2(i)=4*pi*(rvecN(rcross(i,2), icrossvec(i)).^2) * (rvecN(rcross(i,2)+1, icrossvec(i)) - rvecN(rcross(i,2), icrossvec(i)));
    else if rcross(i,1)~=1
        V2(i)=4*pi*(rvecN(rcross(i,2), icrossvec(i)).^2) * (rvecN(rcross(i,2), icrossvec(i)) - rvecN(rcross(i,2)-1, icrossvec(i)));
        end
    end
    
    n1(i)=NumPVir(rcross(i,1))/V1(i);
    n2(i)=NumPVir(rcross(i,2))/V2(i);
    if i~=numel(reac)
        if rcross(i,1)==rcross(i+1,1)
            tdur(i)=tcross(i+1)-tcross(i);
        else
            tdur(i)=0;
        end
    else
        tdur(i)=0;
    end
    FusY(i)=(V1(i)+V2(i))*n1(i)*n2(i)*tdur(i)*reac(i);
end

disp('Total Fusion Yield per cluster')
disp(sum(FusY));

catch
    disp('boo')
end
output.rvecN=rvecN;
output.Q=Q;
output.vvecN=vvecN;
output.NumPVir=NumPVir;
if cross==0
    output.cross=cross;
else
    output.cross=cross;
    output.rcross=rcross;
    output.tcross=tcross;
    output.v1cross=v1cross;
    output.v2cross=v2cross;
end
