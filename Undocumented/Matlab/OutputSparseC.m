function output=OutputSparseC(input)
output=input;
for i=1:numel(input.x)
    output.x{i,1}=sparse(input.x{i,1});
    output.y{i,1}=sparse(input.y{i,1});
    output.z{i,1}=sparse(input.z{i,1});
    output.vx{i,1}=sparse(input.vx{i,1});
    output.vy{i,1}=sparse(input.vy{i,1});
    output.vz{i,1}=sparse(input.vz{i,1});
    output.fEx{i,1}=sparse(input.fEx{i,1});
    output.fEy{i,1}=sparse(input.fEy{i,1});
    output.fEz{i,1}=sparse(input.fEz{i,1});
    output.fBx{i,1}=sparse(input.fBx{i,1});
    output.fBy{i,1}=sparse(input.fBy{i,1});
    output.fBz{i,1}=sparse(input.fBz{i,1});
end