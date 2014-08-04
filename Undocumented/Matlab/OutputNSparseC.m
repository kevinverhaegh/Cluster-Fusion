function output=OutputNSparseC(input)
output=input;
for i=1:numel(output.x)
    output.x{i,1}=full(input.x{i,1});
    output.y{i,1}=full(input.y{i,1});
    output.z{i,1}=full(input.z{i,1});
    output.vx{i,1}=full(input.vx{i,1});
    output.vy{i,1}=full(input.vy{i,1});
    output.vz{i,1}=full(input.vz{i,1});
    output.fEx{i,1}=full(input.fEx{i,1});
    output.fEy{i,1}=full(input.fEy{i,1});
    output.fEz{i,1}=full(input.fEz{i,1});
    output.fBx{i,1}=full(input.fBx{i,1});
    output.fBy{i,1}=full(input.fBy{i,1});
    output.fBz{i,1}=full(input.fBz{i,1});
end