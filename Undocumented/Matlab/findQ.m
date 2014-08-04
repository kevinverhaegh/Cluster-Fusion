function Qi=findQ(Q,rvecN,i,j,k)
%function which returns the value of the enclosed charge of a shell with
%index j at a time with index i of a specie with index k by summing up the
%charges of all enlcosed shells (charge per shell is provided by Q(index
%nr shell, specie).

%essentially, there is a differential equation involving Qi, 
%which means that Qi at time i is actually unknown - and therefore the
%charge should be evaluated of the rvecN of one step earlier (i-1)
if i>1
    i=i-1; 
end

%reshapes the rvecN matrix into a speciexmesh nr matrix containing the
%position of every discretised shell of every specie at time index i
rvecPool=reshape(rvecN(:,i,:),numel(rvecN(:,1,1)),numel(rvecN(1,1,:)));
%retrieves current position of shell index j at time index i of specie index k
r=rvecN(j,i,k); 
%this is a conditional expression, which returns a 1 for every index where
%r>rvecpool and a 0 for every index where r<rvecpool. Hence it returns a 1
%for every enclosed shell (e.q. if the radius of one shell is larger than
%the radius of another shell - it means that this shell encloses the other
%shell).
Bool = r*ones(size(rvecPool))>rvecPool;
Qi=0;
%sum up the charge of every shell which is enclosed of every specie
for jj=1:numel(Bool(:,1))
    for kk=1:numel(Bool(1,:))
       if Bool(jj,kk)==1
           Qi=Qi+Q(kk,jj);
       end
    end
end
%return Qi - the summed up charges

    
