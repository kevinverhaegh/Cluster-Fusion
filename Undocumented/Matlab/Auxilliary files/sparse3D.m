function sparseMat = sparse3D(nonsparse)

[N, index]=min(size(nonsparse));

for i=1:N
    if index==1
        nonsparse(i,:,:)=sparse(nonsparse(i,:,:));
    elseif index==2
        nonsparse(:,i,:)=sparse(nonsparse(:,i,:));
    else
        nonsparse(:,:,i)=sparse(nonsparse(:,:,i));
    end
end

sparseMat=nonsparse;
    