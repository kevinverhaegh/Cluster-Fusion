for i=1:numel(s0)-1
    for j=i+1:numel(s0)
        icross = intersection(rvecN(i,:), rvecN(j,:));
        if icross~=0
            cross=cross+1;
            rcross(cross)=[i, j];
            tcross(cross)=findv(t,icross);
            v1cross(cross)=findv(vvecN(:,i), cross);
            v2cross(cross)=findv(vvecN(:,j),cross);
        end
    end
end