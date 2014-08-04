function AnimateTot(info, tmin, tmax)
 window = 1e-6;
 writerObj = VideoWriter('Animation.avi');
 open(writerObj);
 set(gca,'nextplot','replacechildren');
 set(gcf,'Renderer','zbuffer');
 figure('Units','normalized','Position',[0 0 1 1])
 subplot(2,2,1)
 t = zeros(numel(info),1);
 for i=1:1:numel(info)
     t(i) = info{i,1}.param.time;
 end
 [~,imin] = min(abs(tmin - t));
 [~,imax] = min(abs(tmax - t));
 for i=imin:100:imax
     subplot(2,2,1)
    gscatter(info{i,1}.data.z,info{i,1}.data.y,info{i,1}.data.q,'br','xo')
    axis([-window, window, -window, window])
    title(info{i,1}.param.time/1e-15)
    xlabel('z (m)')
    ylabel('y (m)')
    subplot(2,2,2)
    gscatter(info{i,1}.data.x,info{i,1}.data.z,info{i,1}.data.q,'br','xo')
    axis([-window, window, -window, window])
    title(info{i,1}.param.time/1e-15)
    xlabel('x (m)')
    ylabel('z (m)')
    subplot(2,2,3)
    gscatter(info{i,1}.data.x,info{i,1}.data.y,info{i,1}.data.q,'br','xo')
    axis([-window, window, -window, window])
    title(info{i,1}.param.time/1e-15)
    xlabel('x (m)')
    ylabel('y (m)')
    subplot(2,2,4)
    scatter3(info{i,1}.data.z,info{i,1}.data.y,info{i,1}.data.x, ((info{i,1}.data.q/1.602e-19)+1.1)*5+10, ((info{i,1}.data.q/1.602e-19)+3), 'filled')
    axis([-window, window, -window, window, -window, window])
    title(info{i,1}.param.time/1e-15)
    xlabel('z (m)')
    ylabel('y (m)')
    zlabel('x (m)')
    frame = getframe;
    if sum(abs(info{i,1}.data.x)>window)~=0
    %    window=window*2;
    end
    writeVideo(writerObj,frame);
 end
 close(writerObj);