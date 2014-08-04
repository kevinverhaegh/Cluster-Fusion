function showoutput(input, tb, te)
[~,ib] = min(abs(tb-input.t));
[~,ie] = min(abs(te-input.t));
figure
subplot(2,2,1)
 plot(input.t(ib:ie),input.xge(ib:ie)-2*input.sxge(ib:ie),input.t(ib:ie),input.xge(ib:ie),input.t(ib:ie),input.xge(ib:ie)+2*input.sxge(ib:ie))
 xlabel('t')
 ylabel('Xgem')
 
 subplot(2,2,2)
 plot(input.t(ib:ie),input.yge(ib:ie)-2*input.syge(ib:ie),input.t(ib:ie),input.yge(ib:ie),input.t(ib:ie),input.yge(ib:ie)+2*input.syge(ib:ie))
 xlabel('t')
 ylabel('Ygem')
 
 subplot(2,2,3)
 plot(input.t(ib:ie),input.zge(ib:ie)-2*input.szge(ib:ie),input.t(ib:ie),input.zge(ib:ie),input.t(ib:ie),input.zge(ib:ie)+2*input.szge(ib:ie))
 xlabel('t')
 ylabel('Zgem')
 
 subplot(2,2,4)
 plot3(input.xge(ib:ie),input.yge(ib:ie),input.zge(ib:ie))
 xlabel('Xgem')
 ylabel('Ygem')
 zlabel('Zgem')
 
 figure
 subplot(2,2,1)
 plot(input.t(ib:ie),input.vxge(ib:ie)-2*input.vsxge(ib:ie),input.t(ib:ie),input.vxge(ib:ie),input.t(ib:ie),input.vxge(ib:ie)+2*input.vsxge(ib:ie))
 xlabel('t')
 ylabel('vXgem')
 
 subplot(2,2,2)
 plot(input.t(ib:ie),input.vyge(ib:ie)-2*input.vsyge(ib:ie),input.t(ib:ie),input.vyge(ib:ie),input.t(ib:ie),input.vyge(ib:ie)+2*input.vsyge(ib:ie))
 xlabel('t')
 ylabel('vYgem')
 
 subplot(2,2,3)
 plot(input.t(ib:ie),input.vzge(ib:ie)-2*input.vszge(ib:ie),input.t(ib:ie),input.vzge(ib:ie),input.t(ib:ie),input.vzge(ib:ie)+2*input.vszge(ib:ie))
 xlabel('t')
 ylabel('vZgem')
 
 subplot(2,2,4)
 plot3(input.vxge(ib:ie),input.vyge(ib:ie),input.vzge(ib:ie))
 xlabel('Xgem')
 ylabel('Ygem')
 zlabel('Zgem')