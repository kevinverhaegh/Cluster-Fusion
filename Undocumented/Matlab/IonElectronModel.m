function dy = IonElectronModel(t, y, env, vz0, vx0, eta0, N)
dy = zeros(4,1);
dy(1) = y(2);
dy(2) = (env.omega - env.omega/PhysConst.c * y(4))*dpdeta(1, env.omega*(t - y(3)/PhysConst.c), env, eta0, vz0, vx0) - ...
    ((PhysConst.e^2 * N)/(4*pi*PhysConst.epsilon0))*(y(1))/(sqrt(y(1)^2 + y(3)^2));
dy(3) = y(4);
dy(4) = (env.omega - env.omega/PhysConst.c * y(4))*dpdeta(2, env.omega*(t - y(3)/PhysConst.c), env, eta0, vz0, vx0) - ...
    ((PhysConst.e^2 * N)/(4*pi*PhysConst.epsilon0))*(y(3))/(sqrt(y(1)^2 + y(3)^2));