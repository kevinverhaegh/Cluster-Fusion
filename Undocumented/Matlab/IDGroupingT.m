function output=IDGroupingT(input)

%disp('Grouping different mass/charge combos......')

AMUZcombo=zeros(numel(input)*numel(input{end,1}.data.ID),2);
AMUZcombo(:,1)=input{1,1}.data.m(1);
AMUZcombo(:,2)=input{1,1}.data.q(1);

k=1;
for i=1:numel(input)
    for j=1:numel(input{i,1}.data.ID)
        AMUZcombo(k,1)=input{i,1}.data.m(j);
        AMUZcombo(k,2)=input{i,1}.data.q(j);
        k=k+1;
    end
    %disp(strcat(num2str(100*(i/numel(input))), ' % done'));
end

AMUZcombo=unique(AMUZcombo,'rows');

AMUZcombo(:,1)=AMUZcombo(:,1)/PhysConst.AMU;
AMUZcombo(:,2)=AMUZcombo(:,2)/PhysConst.e;

%disp(strcat(num2str(numel(AMUZcombo(:,1))), ' species have been found:'))

%disp('AMU / Z:')

for i=1:numel(AMUZcombo(:,1))
    %disp(AMUZcombo(i,:))
end
       
%disp('Cathegorising mass/species combos...:')

GroupsID = zeros(numel(input), numel(input{end,1}.data.ID));

for i=1:numel(input)
    for j=1:numel(AMUZcombo(:,1))
        GroupsID(i,1:numel(input{i,1}.data.m))=GroupsID(i,1:numel(input{i,1}.data.m))+j*all(bsxfun(@eq, [input{i,1}.data.m/PhysConst.AMU, input{i,1}.data.q/PhysConst.e], AMUZcombo(j,:)),2)';
    end
   %disp(strcat(num2str(100*(i/numel(input))), ' % done'));
end

%disp('Cathegorisation complete');

%disp('Sorting trajectories based on cathegorisation...');

size=numel(input);
N=numel(input{end,1}.data.ID);
groups = numel(AMUZcombo(:,1));

x = cell(groups,1);
y = cell(groups,1);
z = cell(groups,1);
vx = cell(groups,1);
vy = cell(groups,1);
vz = cell(groups,1);
fEx = cell(groups,1);
fEy = cell(groups,1);
fEz = cell(groups,1);
fBx = cell(groups,1);
fBy = cell(groups,1);
fBz = cell(groups,1);

%disp('Preparing structure...');

for i=1:groups
    x{i,1}=zeros(size,N);
    y{i,1}=zeros(size,N);
    z{i,1}=zeros(size,N);
    vx{i,1}=zeros(size,N);
    vy{i,1}=zeros(size,N);
    vz{i,1}=zeros(size,N);
    fEx{i,1}=zeros(size,N);
    fEy{i,1}=zeros(size,N);
    fEz{i,1}=zeros(size,N);
    fBx{i,1}=zeros(size,N);
    fBy{i,1}=zeros(size,N);
    fBz{i,1}=zeros(size,N);
end

NumberInGroup = zeros(size, groups);
t = zeros(size,1);

%disp('Sorting... This may take a while');

for i=1:size
    t(i) = input{i,1}.param.time;
    zpad=zeros(N-numel(input{i,1}.data.x),1);
    HGroup=GroupsID(i,:);
    Hinput=input{i,1}.data;
    A=numel(input{i,1}.data.x);
    for j=1:numel(AMUZcombo(:,1))
        dum=(HGroup==j)';
        NumberInGroup(i,j)=sum(dum);
        if NumberInGroup(i,j)~=0
        x{j,1}(i,:)=[Hinput.x.*dum(1:A); zpad];
        y{j,1}(i,:)=[Hinput.y.*dum(1:A); zpad];
        z{j,1}(i,:)=[Hinput.z.*dum(1:A); zpad];
        vx{j,1}(i,:)=[Hinput.Bx.*dum(1:A); zpad];
        vy{j,1}(i,:)=[Hinput.By.*dum(1:A); zpad];
        vz{j,1}(i,:)=[Hinput.Bz.*dum(1:A); zpad];
        fEx{j,1}(i,:)=[Hinput.fEx.*dum(1:A); zpad];
        fEy{j,1}(i,:)=[Hinput.fEy.*dum(1:A); zpad];
        fEz{j,1}(i,:)=[Hinput.fEz.*dum(1:A); zpad];
        fBx{j,1}(i,:)=[Hinput.fBx.*dum(1:A); zpad];
        fBy{j,1}(i,:)=[Hinput.fBy.*dum(1:A); zpad];
        fBz{j,1}(i,:)=[Hinput.fBz.*dum(1:A); zpad];
        end
    end
   %disp(strcat(num2str(100*(i/numel(input))), ' % done'));
end

for i=1:groups
    x{i,1}=sparse(x{i,1});
    y{i,1}=sparse(y{i,1});
    z{i,1}=sparse(z{i,1});
    vx{i,1}=sparse(vx{i,1});
    vy{i,1}=sparse(vy{i,1});
    vz{i,1}=sparse(vz{i,1});
    fEx{i,1}=sparse(fEx{i,1});
    fEy{i,1}=sparse(fEy{i,1});
    fEz{i,1}=sparse(fEz{i,1});
    fBx{i,1}=sparse(fBx{i,1});
    fBy{i,1}=sparse(fBy{i,1});
    fBz{i,1}=sparse(fBz{i,1});
end


%disp('Trajectories determined, sorted on cathegorisation.');

%disp('Saving output structure');
            
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
output.GroupsClass=GroupsID;
output.NumberInGroup = NumberInGroup;
output.AMUZcombo = AMUZcombo;
output.param = input{1,1}.param;
