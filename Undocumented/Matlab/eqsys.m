function dydt = eqsys(t,y)
load ionisercalc.mat
%find tcoor
[~,indt]=min(abs(t-tvec));
[~,indr] = min(abs(y(1)-rvec));
Z=Ztab(indt,indr);
[~,indz]=min(abs(Zb-Ztab(indt, indr)));
alpha=alphatab(indr, indz);
    Q=((N*Z)/4)*(2-3*cos(alpha) + cos(alpha)^3);
    dydt = [y(2); ((Q*Z)/y(1)^2)];
    end