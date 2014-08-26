function Energy(input, ibeg)
EkinIons = zeros(numel(input.t),numel(input.AMUZcombo(:,1))-2);
EkinElectrons = zeros(numel(input.t),1);
Epot=zeros(numel(input.t),1);

errors=find(sum((input.x{1,1}~=0)')'~=input.NumberInGroup(:,1));
try
for i=ibeg:numel(input.t)
    if input.NumberInGroup(i,1)~=0 && sum(i==errors)~=1
    EkinElectrons(i)=sum(1./sqrt(1-(input.vx{1,1}(i,:).^2 + input.vy{1,1}(i,:).^2 + input.vz{1,1}(i,:).^2)))*PhysConst.me*PhysConst.c^2;
    for j=3:numel(input.AMUZcombo(:,1))
        EkinIons(i,j-2)=sum(1./sqrt(1-(input.vx{j,1}(i,:).^2 + input.vy{j,1}(i,:).^2 + input.vz{j,1}(i,:).^2)))*input.AMUZcombo(j,1)*PhysConst.AMU*PhysConst.c^2;
    end
    %group all info of charged particles for potential energy calculation
    partx=[];
    party=[];
    partz=[];
    partq=[];
    for j=1:numel(input.AMUZcombo(:,1))
        if j~=2
        dum=[partx, full(input.x{j,1}(i,input.GroupsClass(i,:)==j))];
        partx=dum;
        dum=[party, full(input.y{j,1}(i,input.GroupsClass(i,:)==j))];
        party=dum;
        dum=[partz, full(input.z{j,1}(i,input.GroupsClass(i,:)==j))];
        partz=dum;
        dum=[partq, input.AMUZcombo(j,2)*ones(1,input.NumberInGroup(i,j))];
        partq=dum;
        end
    end

    Epoth=0;
    for j=1:numel(partx)-1
        Epoth=Epoth + ((PhysConst.e^2)/(4*pi*PhysConst.epsilon0))*sum((partq(j).*partq(j+1:end))./sqrt((partx(j)-partx(j+1:end)).^2 + (party(j)-party(j+1:end)).^2 + (partz(j)-partz(j+1:end)).^2));
    end
    Epot(i)=Epoth;
    end
    if sum(i==errors)==1
        EkinElectrons(i)=EkinElectrons(i-1);
        EkinIons(i,:)=EkinIons(i-1,:);
        Epot(i,:)=Epot(i-1,:);
    end
end
catch
    disp('boo')
disp(i)
end

if numel(input.AMUZcombo(:,1))>3
Etot = EkinElectrons + Epot + sum(EkinIons')';
else
    Etot = EkinElectrons + Epot + EkinIons;
end
plot(input.t, EkinIons, input.t, EkinElectrons, input.t, Epot, input.t, Etot)
