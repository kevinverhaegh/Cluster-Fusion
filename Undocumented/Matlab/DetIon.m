function t=DetIon(G,info)
for i=1:numel(info)
  t=info{i,1}.param.time;
    if numel(info{i,1}.data.x)==G
        break;
    end
end