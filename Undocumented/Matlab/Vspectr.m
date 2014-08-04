function Vspectru = Vspectr(tind,vvecN,NumPVir)
size=cumsum(floor(NumPVir/100));
Vspectru=zeros(1,size(end));
for j=2:numel(NumPVir)
   if size(j-1)~=0
   Vspectru(size(j-1):size(j))=ones(1,numel(size(j-1):size(j))) .*vvecN(j-1,tind);
   end
end

        
        