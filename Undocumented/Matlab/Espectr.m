function Vspectr(t,vvecN,NumPVir)
size=cumsum(floor(NumPVir));
Espectr=zeros(numel(t),size(end));
for i=1:numel(t)
    for j=2:numel(NumPVir)
        if size(j-1)~=0
        Espectr(i,size(j-1):size(j))=ones(1,size(j-1):size(j)) .*vvecN(j,i);
        end
    end
end
        
        