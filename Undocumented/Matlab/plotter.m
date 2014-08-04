close all
subplot(3,1,1)
plot(output.t(45000:75000)/1e-15,1e8 *full(output.x{1,1}(45000:75000,101:200)))
axis([90 150 -20 20])
xlabel('t (fs)', 'FontName','Helvetica')
ylabel('x (10^{-8} m)', 'FontName','Helvetica')
subplot(3,1,2)
plot(output.t(45000:75000)/1e-15,1e8 *full(output.y{1,1}(45000:75000,101:200)))
axis([90 150 -6 6])
xlabel('t (fs)', 'FontName','Helvetica')
ylabel('y (10^{-8} m)', 'FontName','Helvetica')
subplot(3,1,3)
plot(output.t(45000:75000)/1e-15,1e8 *full(output.z{1,1}(45000:75000,101:200)))
axis([90 150 -6 6])
xlabel('t (fs)', 'FontName','Helvetica')
ylabel('z (10^{-8} m)', 'FontName','Helvetica')
ax=samexaxis('xmt','on','ytac','join','yld',1);
set( ax(3)                       , ...
    'FontName'   , 'Helvetica' );
title(ax(1),'Electron Trajectories','FontName','Helvetica','FontSize',12,'FontWeight','bold')
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 'C:\Users\Kevin Verhaegh\Dropbox\Documents\Graduate project Cluster Fusion\Report\figures\C3Coul_ElectronTrajectories.eps'

figure
subplot(3,1,1)
plot(output.t/1e-15,1e8 * full(output.x{3,1}(:,1:100)))
xlabel('t (fs)')
ylabel('x (10^{-8} m)')
axis([0 1000 -5 5])
subplot(3,1,2)
plot(output.t/1e-15,1e8 * full(output.y{3,1}(:,1:100)))
xlabel('t (fs)')
ylabel('y (10^{-8} m)')
axis([0 1000 -5 5])
subplot(3,1,3)
plot(output.t/1e-15,1e8 * full(output.z{3,1}(:,1:100)))
xlabel('t (fs)')
ylabel('z (10^{-8} m)')
axis([0 1000 -5 5])
ax=samexaxis('xmt','on','ytac','join','yld',1);
title(ax(1),'Ion Trajectories','FontName','Helvetica','FontSize',12,'FontWeight','bold')
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 'C:\Users\Kevin Verhaegh\Dropbox\Documents\Graduate project Cluster Fusion\Report\figures\C3Coul_IonTrajectories.eps'
close;

figure
subplot(3,1,1)
plot(output.t(45000:75000)/1e-15,1e-6 * PhysConst.c*full(output.vx{1,1}(45000:75000,101:200)))
axis([90 150 -12.5 12.5])
xlabel('t (fs)')
ylabel('v_x (10^6 m/s)')
subplot(3,1,2)
plot(output.t(45000:75000)/1e-15,1e-6 * PhysConst.c*full(output.vy{1,1}(45000:75000,101:200)))
axis([90 150 -4 4])
xlabel('t (fs)')
ylabel('v_y (10^6 m/s)')
subplot(3,1,3)
plot(output.t(45000:75000)/1e-15,1e-6 * PhysConst.c*full(output.vz{1,1}(45000:75000,101:200)))
axis([90 150 -4 4])
xlabel('t (fs)')
ylabel('v_z (10^6 m/s)')
samexaxis('xmt','on','ytac','join','yld',1);
ax=samexaxis('xmt','on','ytac','join','yld',1);
title(ax(1),'Electron Velocity Traces','FontName','Helvetica','FontSize',12,'FontWeight','bold')
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 'C:\Users\Kevin Verhaegh\Dropbox\Documents\Graduate project Cluster Fusion\Report\figures\C3Coul_ElectronVelocities.eps'
close;

figure
subplot(3,1,1)
plot(output.t / 1e-15,1e-4 * PhysConst.c*full(output.vx{3,1}(:,1:100)))
axis([0 1000 -6 6])
xlabel('t (fs)')
ylabel('v_x (10^4 m/s)')
subplot(3,1,2)
plot(output.t / 1e-15,1e-4 * PhysConst.c*full(output.vy{3,1}(:,1:100)))
xlabel('t (fs)')
ylabel('v_y (10^4 m/s)')
axis([0 1000 -6 6])
subplot(3,1,3)
plot(output.t / 1e-15,1e-4 * PhysConst.c*full(output.vz{3,1}(:,1:100)))
xlabel('t (fs)')
ylabel('v_z (10^4 m/s)','FontName','Helvetica')
axis([0 1000 -6 6])
ax=samexaxis('xmt','on','ytac','join','yld',1);
title(ax(1),'Ion Velocity Traces','FontName','Helvetica','FontSize',12,'FontWeight','bold')
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 'C:\Users\Kevin Verhaegh\Dropbox\Documents\Graduate project Cluster Fusion\Report\figures\C3Coul_IonVelocities.eps'
close;

