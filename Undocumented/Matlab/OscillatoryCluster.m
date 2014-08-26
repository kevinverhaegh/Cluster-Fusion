function xx=OscillatoryCluster(E,rho,x0,R)

omegaL=2.354564e+15;

A = E*PhysConst.e/PhysConst.me;
B = rho*PhysConst.e^2/(3*PhysConst.epsilon0*PhysConst.me);

% Phase portrait for the Duffing system.
deq=@(t,x) [x(2);-B*CoulFC(x(1),R)+A*sin(omegaL*t)]; 
options=odeset('RelTol',1e-5);
[t,xx]=ode45(deq,0:(2*pi / omegaL)*0.01:5e-14,[x0,0],options);