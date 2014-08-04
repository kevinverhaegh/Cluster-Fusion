function output=ShockShellModelGeneral(rho, R0, Z, M, t, r0t, ReacCross)

%input: for each specie rho(r,specie), R0(specie), Z(specie), M(specie), 
%t(time), r0t(r), ReacCross(specie, specie). ReacCross is a square matrix
%indicating which species can fuse with each other, and which type of
%fusion reaction this is (DD, DT, etc.).

sizevec=size(rho);
NrSpecies=min(sizevec);

disp('Number of species:')
disp(NrSpecies)

N=zeros(NrSpecies,1);

for i=1:NrSpecies
N(i)=trapz(r0t, rho(:,i).*4.*pi.*r0t.^2);
end

disp('Amount of ions:')
disp(N)

densavg=sum(Z.*N)./(4/3 * pi .* R0.^3); 

t0 = sqrt(3.*M.*PhysConst.epsilon0./(densavg .* Z * PhysConst.e^2));

disp('CE time scales:')
disp(t0)

s0=r0t./R0(i);
tau=zeros(NrSpecies,numel(t));

for i=1:NrSpecies
tau(i,:) = t./t0(i);
end

Q=zeros(NrSpecies,numel(r0t));
for i=1:NrSpecies
Q(i,:)=[0 diff(cumtrapz(r0t, 4.*pi.*Z(i).*rho(:,i).*r0t.^2))'];
end
rvecN=zeros(numel(rho(:,1)), numel(t),NrSpecies);
Q=Q./(4/3 .* pi .* R0(1).^3 .* densavg(1));
for i=1:NrSpecies
rvecN(:,1,i)=r0t;
end

%try
%fun=@(x,Qi,s0i,taui) (sqrt(x*(x-1)) + log(sqrt(x) + sqrt(x-1)) - taui*sqrt(2*Qi/(s0i^3)));
fun=@(x,Qi,s0i,taui,a) (sqrt(x*(x-a)) + 0.5*a*log(2*sqrt(x*(x-a)) + x - a) - taui*sqrt(2*Qi/(s0i^3)));
for i=1:numel(t)
    for k=1:NrSpecies
        taui=tau(k,i);
        for j=2:numel(s0)
            Qi=findQ(Q,rvecN,i,j,k);
            s0i=s0(j);
            if j>1
            a=Qi/sum(sum(Q(:,1:j)));
            else
                a=1;
            end
            %specfun = @(x) fun(x, Qi, s0i, taui) + (x>1);
            specfun = @(x) fun(x, Qi, s0i, taui, a) + (x>a);
            opts = optimset('Diagnostics','off', 'Display','off', 'TolFun', 1e-9, 'TolX', 1e-9);
            rvecN(j,i,k)=r0t(j)*real(fsolve(specfun, a+1e-6,opts));
        end
    end
end

vvecN=zeros(numel(rvecN(:,1,1)),numel(rvecN(1,:,1))-1,numel(rvecN(1,1,:)));

for i=1:NrSpecies
figure
plot(t,rvecN(:,:,i)')
end
for i=1:NrSpecies
vvecN(:,:,i)=diff(rvecN(:,:,i)')'/(t(2)-t(1));
end
%each virtual particle attributes to a shell, which represents a number of
%particles, calculated here

%WORKS UNTIL HERE !!!!!!!!!!!!!!!!!!!
NumPVir=zeros(numel(r0t),NrSpecies);
for i=1:NrSpecies
NumPVirH=cumtrapz([0 r0t], 4*pi*[rho(1,i) rho(:,i)].*[0 r0t].^2);
NumPVir(:,i)=[NumPVirH(1) diff(NumPVirH)];
end

cross=0;
totcross=0;
%search for collisions, a collision is an intersection of the various lines
%of various radia, first search which collisions can occur.

ColVec=zeros(NrSpecies^2);
ColMat=zeros(NrSpecies^2,2);
Col=0;
%First count how many species can collide, don't count doubles
for k=1:NrSpecies
    for l=1:NrSpecies-k+1
        if ReacCross(k,l)~=0
            Col=Col+1;
            ColVec(Col)=ReacCross(k,l);
            ColMat(Col,:)=[k,l]
        end
    end
end

%Count collisions for each possible collisionable specie

icrossvec=zeros(numel(r0t).^2,Col);
tcross=0.*icrossvec;
v1cross=0.*icrossvec;
v2cross=0.*icrossvec;
rcross=zeros(numel(r0t).^2,2,Col);
reac=0.*v1cross;
Erel=0.*v1cross;

for k=1:Col
    for i=1:numel(s0)-1
        for j=i+1:numel(s0)
            icross = intersection(rvecN(i,:,ColMat(k,1)), rvecN(j,:,ColMat(k,2)));
            if icross~=0
                cross=cross+1;
                icrossvec(cross,k)=round(icross);
                rcross(cross, :,k)=[i, j];
                tcross(cross,k)=findv(t,icross);
                v1cross(cross,k)=findv(vvecN(:,i,ColMat(k,1)), icross);
                v2cross(cross,k)=findv(vvecN(:,j,ColMat(k,2)),icross);
            end
        end
    end

    %calculate fusion reactions
    %load DD crossections
    if ColVec(k)==1
    load('DD_crossection.mat');
    else if ColVec(k)==2
            load('DT_crossection.mat');
        end
    end
    Erel(:,k) = 0.5 * M * abs(v1cross(:,k) - v2cross(:,k)).^2 ./ PhysConst.e;
    for i=1:numel(reac(:,1))
        if Erel(i,k)<20
        reac(i,k)=0;
        else
        [~,r]=min(abs(Erel(i,k)-EnDDsp));
        reac(i,k)=DDsp(r)*(10^-28)*abs(v1cross(i,k)-v2cross(i,k));
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
    V1(i)=abs(4*pi*(rvecN(rcross(i,1), icrossvec(i)).^2) * (rvecN(rcross(i,1)+1, icrossvec(i)) - rvecN(rcross(i,1), icrossvec(i))));
    else if rcross(i,1)~=1
        V1(i)=abs(4*pi*(rvecN(rcross(i,1), icrossvec(i)).^2) * (rvecN(rcross(i,1), icrossvec(i)) - rvecN(rcross(i,1)-1, icrossvec(i))));
        end
    end
    
    if rcross(i,1)~=numel(rvecN(:,2)) && rcross(i,2)~=numel(rvecN(:,1))
    V2(i)=abs(4*pi*(rvecN(rcross(i,2), icrossvec(i)).^2) * (rvecN(rcross(i,2)+1, icrossvec(i)) - rvecN(rcross(i,2), icrossvec(i))));
    else if rcross(i,1)~=1
        V2(i)=abs(4*pi*(rvecN(rcross(i,2), icrossvec(i)).^2) * (rvecN(rcross(i,2), icrossvec(i)) - rvecN(rcross(i,2)-1, icrossvec(i))));
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
end
