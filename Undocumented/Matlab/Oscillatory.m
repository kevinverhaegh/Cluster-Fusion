function xx=Oscillatory(E,D)

Z=1;
omegaL=2.354564e+15;

A = E*PhysConst.e/PhysConst.me;
B = Z*PhysConst.e^2/(4*pi*PhysConst.epsilon0*PhysConst.me);
C=2e-10;

% Phase portrait for the Duffing system.
deq=@(t,x) [x(2);-B*CoulF(x(1),C)+A*cos(omegaL*t)]; 
options=odeset('RelTol',1e-5);
[t,xx]=ode45(deq,0:(2*pi / omegaL)*0.01:5e-14,[D,0],options);

%subplot(2,1,1)
plot(xx(:,1),xx(:,2))
%subplot(2,1,2)
%plot(t,xx(:,1))

%fsize=15;
axis([-1e-9 1e-9 -2.5e6 2.5e6])
%xlabel('x','FontSize',fsize)
%ylabel('y','FontSize',fsize)

% End of Programs_14c.