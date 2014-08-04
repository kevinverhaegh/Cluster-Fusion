    tsize=3000;
    x = zeros(tsize,1,2,100);
    y = zeros(tsize,1,2,100);
    z = zeros(tsize,1,2,100);
    vx = zeros(tsize,1,2,100);
    vy = zeros(tsize,1,2,100);
    vz = zeros(tsize,1,2,100);
parfor i=1:100
    if i< 10
    info=InputFix(loadgdf(strcat('UCP00',num2str(i), '.gdf')));
    else if i < 100
            info=InputFix(loadgdf(strcat('UCP0',num2str(i), '.gdf')));
        end
    end
    if i==100
        info=InputFix(loadgdf(strcat('UCP',num2str(i),'.gdf')));
    end
    
    [x(:,:,:,i), y(:,:,:,i), z(:,:,:,i), vx(:,:,:,i), vy(:,:,:,i), vz(:,:,:,i), t, groupsID] = Trajectories(info);
end