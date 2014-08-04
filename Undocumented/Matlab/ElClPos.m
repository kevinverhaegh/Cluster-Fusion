function [Te,output] = ElClPos(input)
t=zeros(numel(input),1);
xge=zeros(numel(input),1);
yge=zeros(numel(input),1);
zge=zeros(numel(input),1);
sxge = zeros(numel(input),1);
syge=zeros(numel(input),1);
szge=zeros(numel(input),1);
xgi=zeros(numel(input),1);
ygi=zeros(numel(input),1);
zgi=zeros(numel(input),1);
sxgi = zeros(numel(input),1);
sygi=zeros(numel(input),1);
szgi=zeros(numel(input),1);
Te = cell(numel(input),1);
Teg=zeros(numel(input),1);
IonComp = 0;
EleComp = 0;

 for i=1:numel(input)
     if IonComp ~= 1
         IonCount = 0;
         for j=1:numel(input{i,1}.data.q)
             if input{i,1}.data.q(j) > 0
                 IonCount = IonCount + 1;
             end
         end
         if IonCount == input{1,1}.param.nps
             IonComp = 1;
         end
     end
     if EleComp ~= 1
         EleCount = 0;
         for j=1:numel(input{i,1}.data.q)
             if input{i,1}.data.q(j) < 0
                 EleCount = EleCount + 1;
             end
         end
         if EleCount == input{1,1}.param.nps
             EleComp = 1;
             tioncomp = input{i,1}.param.time;
         end
     end
     
     t(i) = input{i,1}.param.time;
     
     xge(i) = mean(input{i,1}.data.x(1+IonCount:IonCount+EleCount));
     xgi(i) = mean(input{i,1}.data.x(1:IonCount));
     yge(i) = mean(input{i,1}.data.y(1+IonCount:IonCount+EleCount));
     ygi(i) = mean(input{i,1}.data.y(1:IonCount));
     zge(i) = mean(input{i,1}.data.z(1+IonCount:IonCount+EleCount));
     zgi(i) = mean(input{i,1}.data.z(1:IonCount));

     sxge(i) = std(input{i,1}.data.x(1+IonCount:IonCount+EleCount));
     sxgi(i) = std(input{i,1}.data.x(1:IonCount));
     syge(i) = std(input{i,1}.data.y(1+IonCount:IonCount+EleCount));
     sygi(i) = std(input{i,1}.data.y(1:IonCount));
     szge(i) = std(input{i,1}.data.z(1+IonCount:IonCount+EleCount));
     szgi(i) = std(input{i,1}.data.z(1:IonCount));
     
     Te{i,1} = TeDistr(input{i,1}.data.Bx(1+IonCount:IonCount+EleCount),input{i,1}.data.By(1+IonCount:IonCount+EleCount),input{i,1}.data.Bz(1+IonCount:IonCount+EleCount));
     Teg(i) = mean(Te{i,1});
     disp((i/numel(input)) * 100)
 end
 
 subplot(2,2,1)
 plot(t,xge-2*sxge,t,xge,t,xge+2*sxge)
 xlabel('t')
 ylabel('Xgem')
 
 subplot(2,2,2)
 plot(t,yge-2*syge,t,yge,t,yge+2*syge)
 xlabel('t')
 ylabel('Ygem')
 
 subplot(2,2,3)
  plot(t,zge-2*szge,t,zge,t,zge+2*szge)
 xlabel('t')
 ylabel('Xgem')
 
 subplot(2,2,4)
 plot3(xge,yge,zge)
 xlabel('Xgem')
 ylabel('Ygem')
 zlabel('Zgem')
 
 output = struct('Teg', Teg,'tioncomp', tioncomp, 't', t,'xge',xge,'xgi',xgi, 'yge',yge,'ygi',ygi,'zge',zge,'zgi',zgi,'sxge',sxge,'sxgi',sxgi, 'syge',syge,'sygi',sygi,'szge',szge,'szgi',szgi);