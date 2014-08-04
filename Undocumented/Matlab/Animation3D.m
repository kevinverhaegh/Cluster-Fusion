function M=Animation3D(info)
close all
for i=1:numel(info)
    scatter3(info{i,1}.data.z,info{i,1}.data.y,info{i,1}.data.x, ((info{i,1}.data.q/1.602e-19)+1.1)*5+10, ((info{i,1}.data.q/1.602e-19)+3), 'filled')
    axis([-2e-7, 2e-7, -2e-7, 2e-7, -2e-7, 2e-7])
    title(info{i,1}.param.time/1e-15)
    xlabel('z (m)')
    ylabel('y (m)')
    zlabel('x (m)')
    M(i) = getframe();
end