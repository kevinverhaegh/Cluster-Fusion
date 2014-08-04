function rpred=CoulombExpModel(input, specie)
rs=max(sqrt(input.x{specie,1}.^2 + input.y{specie,1}.^2 + input.z{specie,1}.^2)');
dum=max(sqrt(input.x{2,1}.^2 + input.y{2,1}.^2 + input.z{2,1}.^2)');
r0=dum(1);
t=input.t;
wp0=sqrt((input.param.dens.*PhysConst.e^2)/(input.AMUZcombo(specie,1).*PhysConst.AMU.*PhysConst.epsilon0));

tau=(t-t(1)).*wp0;
[~,idelay]=min(abs((input.NumberInGroup(:,3)==input.NumberInGroup(1,2)) - 1));
tdelay=t(idelay);
%assume spheriod

rtesta=1.000000001:1e-4:10;
rtestb=11:1:1e3;

rtest=zeros(1,numel(rtesta)+numel(rtestb));

rtest(1:numel(rtesta))=rtesta;
rtest(numel(rtesta)+1:numel(rtesta)+numel(rtestb))=rtestb;


testfunc=cumtrapz(rtest, sqrt(rtest./(rtest-1)));

rpred=0.*tau;

for i=1:numel(tau)
[~,dum]=min(abs(tau(i)-sqrt(3/2).*testfunc));
rpred(i)=rtest(dum)*r0;
end

figure
plot(t,rs,t+tdelay,rpred)
aap=1;