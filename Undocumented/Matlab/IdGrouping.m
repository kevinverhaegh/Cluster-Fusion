function output=IdGrouping(input,Parallel, AMUZcombo)

if Parallel==1
parpool
end

if nargin==2
    
disp('Grouping different mass/charge combos......')

AMUZcombo=zeros(numel(input)*numel(input{end,1}.data.ID),2);
AMUZcombo(:,1)=input{1,1}.data.m(1)/PhysConst.AMU;
AMUZcombo(:,2)=input{1,1}.data.q(1)/PhysConst.e;

k=1;
for i=1:numel(input)
    for j=1:numel(input{i,1}.data.ID)
        AMUZcombo(k,1)=input{i,1}.data.m(j)/PhysConst.AMU;
        AMUZcombo(k,2)=input{i,1}.data.q(j)/PhysConst.e;
        k=k+1;
    end
    disp(strcat(num2str(100*(i/numel(input))), ' % done'));
end

AMUZcombo=unique(AMUZcombo,'rows');

disp(strcat(num2str(numel(AMUZcombo(:,1))), ' species have been found:'))

disp('AMU / Z:')

for i=1:numel(AMUZcombo(:,1))
    disp(AMUZcombo(i,:))
end
            
        
disp('Cathegorising mass/species combos...:')

GroupsID = zeros(numel(input), numel(input{end,1}.data.ID));

for i=1:numel(input)
    for j=1:numel(AMUZcombo(:,1))
        GroupsID(i,1:numel(input{i,1}.data.m))=GroupsID(i,1:numel(input{i,1}.data.m))+j*ismember([input{i,1}.data.m/PhysConst.AMU, input{i,1}.data.q/PhysConst.e], AMUZcombo(j,:), 'rows')';
    end
   disp(strcat(num2str(100*(i/numel(input))), ' % done'));
end

disp('Cathegorisation complete');
else if nargin==2
        GroupsID = zeros(numel(input), numel(input{end,1}.data.ID));
        for i=1:numel(input)
            for j=1:numel(AMUZcombo(:,1))
                 GroupsID(i,1:numel(input{i,1}.data.m))=GroupsID(i,1:numel(input{i,1}.data.m))+j*ismember([input{i,1}.data.m/PhysConst.AMU, input{i,1}.data.q/PhysConst.e], AMUZcombo(j,:), 'rows')';
            end
            disp(strcat(num2str(100*(i/numel(input))), ' % done'));
        end
    end
end

disp('Sorting trajectories based on cathegorisation...');

size=numel(input);
N=numel(input{end,1}.data.ID);
groups = numel(AMUZcombo(:,1));

NumberInGroup = zeros(size, groups);
t = zeros(size,1);

parfor i=1:size
    t(i) = input{i,1}.param.time;
    for j=1:groups
            NumberInGroup(i,j)=sum(ismember(GroupsID(i,:), j));
    end
end

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
IDMat = cell(groups,1);

parfor i=1:groups
    x{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    y{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    z{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    vx{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    vy{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    vz{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    fEx{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    fEy{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    fEz{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    fBx{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    fBy{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    fBz{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
    IDMat{i,1}=sparse(zeros(numel(t),max(NumberInGroup(:,i))));
end

for i=1:size
    for j=1:numel(input{i,1}.data.ID)
        [~,indS]=find([input{i,1}.data.m(j)/PhysConst.AMU input{i,1}.data.q(j)/PhysConst.e] == AMUZcombo);
        [~,ind]=find(IDMat{indS,1}(i,:));
        indP=min(ind);
        x{indS,1}(i,indP)=input{i,1}.data.x(j);
        y{indS,1}(i,indP)=input{i,1}.data.y(j);
        z{indS,1}(i,indP)=input{i,1}.data.z(j);
        vx{indS,1}(i,indP)=input{i,1}.data.vx(j);
        vy{indS,1}(i,indP)=input{i,1}.data.vy(j);
        vz{indS,1}(i,indP)=input{i,1}.data.vz(j);
        fEx{indS,1}(i,indP)=input{i,1}.data.fEx(j);
        fEy{indS,1}(i,indP)=input{i,1}.data.fEy(j);
        fEz{indS,1}(i,indP)=input{i,1}.data.fEz(j);
        fBx{indS,1}(i,indP)=input{i,1}.data.fBx(j);
        fBy{indS,1}(i,indP)=input{i,1}.data.fBy(j);
        fBz{indS,1}(i,indP)=input{i,1}.data.fBz(j);
        IDMat{indS,1}(i,indP)=input{i,1}.data.ID(j);
    end
end



disp('Trajectories determined, sorted on cathegorisation.');

disp('Saving output structure');
            
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

if Parallel == 1
delete(gcp);
end