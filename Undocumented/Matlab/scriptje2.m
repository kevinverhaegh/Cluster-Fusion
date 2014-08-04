FcRet = zeros(9, numel(outputCoul.t), 3);
FcDir = zeros(9, numel(outputCoul.t), 3);
FcRetNRad = zeros(9, numel(outputCoul.t),3);
FcRetTot = zeros(numel(outputCoul.t),3);
FcDirTot = zeros(numel(outputCoul.t),3);
FcRetNRadTot = zeros(numel(outputCoul.t),3);
for i=1:9
outputs=CompareRetard3(trajec(1), trajec(i+1), outputCoul.t);
FcRetTot = FcRetTot + outputs.FcRet;
FcDirTot = FcDirTot + outputs.FcInst;
FcRetNRadTot = FcRetNRadTot + outputs.FcRetNoRad;
FcRet(i,:,:) = outputs.FcRet;
FcDir(i,:,:) = outputs.FcInst;
FcRetNRad(i,:,:) = outputs.FcRetNoRad;
end