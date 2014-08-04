function info=ElectronField(type, x0,y0,z0,env,vx0,vy0,vz0,t)
m=PhysConst.me;
        e=PhysConst.e;
        c=PhysConst.c;
        E0=env.A;
        omega=env.omega;
        eta1=omega*(t(1)-(z0/c));
        %calculate pondermotive quantities
        
        %pondermotive energy
        
        ponden = e^2 * env.A /(4 * m *omega^2);
        pondposmax = e * env.A /(4 * m *omega^2);
        pondermotive = struct('ponden',ponden,'pondposmax',pondposmax);
  
if strcmp(type, 'PlaneTravelPhi')==1
    disp('Using traveling plane wave, travelling in z direction, polarization in y direction')
    if vy0==0
        disp('vy0 = 0, using Peters theory, assuming t=0')
        %assuming steps in t don't vary, assuming t=0
        vy = 0.*t;
        vy(1) = 0;
        vz = 0.*t;
        vz(1) = vz0;
        vx=0.*t;
        vx(1)=vx0;
        gamma = zeros(numel(t) - 1,1);
        gamma(1) = 1/sqrt(1 - ((vz0^2 + vx0^2)/PhysConst.c^2));
        x=0.*t;
        x(1)=x0;
        y=0.*t;
        y(1)=y0;
        z=0.*t;
        z(1)=z0;
        pvz=0.*t;
        pvz(1)=gamma(1)*m*vz0;
        pvx=0.*t;
        pvx(1)=gamma(1)*m*vx0;
        pvzd=0.*t;
        pvzd(1)=gamma(1)*m*vz0;
        pvyd=0.*t;
        pvyd(1)=gamma(1)*m*vy0;
        pvxd=0.*t;
        pvy=0.*t;
        vab=0.*t;
        vab(1) = vz0^2 + vy0^2 /PhysConst.c^2;
        eta=0.*t;
        eta(1)=eta1;
        gammad = zeros(numel(t) - 1,1);
        gammad(1) = 1/sqrt(1 - ((vz0^2 + vy0^2+vx0^2)/PhysConst.c^2));
        xd=0.*t;
        xd(1)=x0;
        yd=0.*t;
        yd(1)=y0;
        zd=0.*t;
        zd(1)=z0;
        vyd = 0.*t;
        vyd(1) = vy0;
        vzd = 0.*t;
        vzd(1) = vz0;
        vxd=0.*t;
        A=-e*E0/(omega);
        phi = env.phi;
        %creating polynomial to find root vz
        for i=2:numel(t)
            %numerical approximation, taking gamma of previous timestep
            %Finding vz using gamma(i-1)
            eta(i)=omega*(t(i)-(z(i-1)/c));
            pvz(i) =  m * gamma(1) * vz(1) + ((0.5*(A^2)*((sin(eta(i)+phi) - sin(phi))^2) + (A*pvx(1)*(sin(eta(i)+phi)-sin(phi)))) / (m*gamma(1)*(c-vz(1))));
            pvx(i) = m*gamma(1)*vx(1) + (A*(sin(eta(i)+phi)-sin(phi)));
            vab(i)=sqrt(((pvz(i)^2 + pvx(i)^2))/(c^2 * m^2 + (pvz(i)^2 + pvx(i)^2))); 
            gamma(i) = 1/(sqrt((1 + (vab(i)))*(1 - (vab(i)))));
            vx(i) = pvx(i)/(m*gamma(i));
            vz(i) = pvz(i)/(m*gamma(i));
            y(i) = 0;
            z(i) = z(i-1) + vz(i)*(t(i) - t(i-1));
            x(i) = x(i-1) + vx(i)*(t(i) - t(i-1));
            vy(i) = 0;
            pvy(i) = 0;
            pvzd(i) = pvz(1) + (0.5/(m*(gamma(1)*(c-vz0))))*(0.5*(A^2 * (1 + 2*(sin(phi)^2))) - 2*pvx(1)*A*sin(phi));
            pvxd(i) = pvx(1) - A *sin(phi);
            pvyd(i) = 0;
            %disp((i/numel(t))*100)
        end
    end
end
    
    
if strcmp(type, 'PlaneTravel')==1
    disp('Using traveling plane wave, travelling in z direction, polarization in y direction')
    if vy0==0
        disp('vy0 = 0, using Peters theory, assuming t=0')
        %assuming steps in t don't vary, assuming t=0
        vy = 0.*t;
        vy(1) = 0;
        vz = 0.*t;
        vz(1) = vz0;
        vx=0.*t;
        vx(1)=vx0;
        gamma = zeros(numel(t) - 1,1);
        gamma(1) = 1/sqrt(1 - ((vz0^2 + vx0^2)/PhysConst.c^2));
        x=0.*t;
        x(1)=x0;
        y=0.*t;
        y(1)=y0;
        z=0.*t;
        z(1)=z0;
        pvz=0.*t;
        pvz(1)=gamma(1)*m*vz0;
        pvx=0.*t;
        pvx(1)=gamma(1)*m*vx0;
        pvzd=0.*t;
        pvzd(1)=gamma(1)*m*vz0;
        pvyd=0.*t;
        pvyd(1)=gamma(1)*m*vy0;
        pvxd=0.*t;
        pvy=0.*t;
        vab=0.*t;
        vab(1) = vz0^2 + vy0^2 /PhysConst.c^2;
        eta=0.*t;
        eta(1)=eta1;
        gammad = zeros(numel(t) - 1,1);
        gammad(1) = 1/sqrt(1 - ((vz0^2 + vy0^2+vx0^2)/PhysConst.c^2));
        xd=0.*t;
        xd(1)=x0;
        yd=0.*t;
        yd(1)=y0;
        zd=0.*t;
        zd(1)=z0;
        vyd = 0.*t;
        vyd(1) = vy0;
        vzd = 0.*t;
        vzd(1) = vz0;
        vxd=0.*t;
        A=-e*E0/(omega);
        %creating polynomial to find root vz
        for i=2:numel(t)
            %numerical approximation, taking gamma of previous timestep
            %Finding vz using gamma(i-1)
            eta(i)=omega*(t(i)-(z(i-1)/c)); 
            pvz(i) =  m * gamma(1) * vz(1) + ((0.5*(A^2)*((sin(eta(i)) - sin(eta(1)))^2) + (A*pvx(1)*(sin(eta(i))-sin(eta(1))))) / (m*gamma(1)*(c-vz(1))));
            pvx(i) = m*gamma(1)*vx(1) + (A*(sin(eta(i))-sin(eta(1))));
            vab(i)=sqrt(((pvz(i)^2 + pvx(i)^2))/(c^2 * m^2 + (pvz(i)^2 + pvx(i)^2))); 
            gamma(i) = 1/(sqrt((1 + (vab(i)))*(1 - (vab(i)))));
            vx(i) = pvx(i)/(m*gamma(i));
            vz(i) = pvz(i)/(m*gamma(i));
            y(i) = 0;
            z(i) = z(i-1) + vz(i)*(t(i) - t(i-1));
            x(i) = x(i-1) + vx(i)*(t(i) - t(i-1));
            vy(i) = 0;
            pvy(i) = 0;
            pvzd(i) = pvz(1) + (0.5/(m*(gamma(1)*(c-vz0))))*(0.5*(A^2 * (1 + 2*(sin(eta(1))^2))) - 2*pvx(1)*A*sin(eta(1)));
            pvxd(i) = pvx(1) - A *sin(eta(1));
            pvyd(i) = 0;
            %disp((i/numel(t))*100)
        end
    end
    
elseif strcmp(type, 'GeneralTravel')==1
        
        vy = 0.*t;
        vy(1) = vy0;
        vz = 0.*t;
        vz(1) = vz0;
        vx=0.*t;
        vx(1)=vx0;
        gamma = zeros(numel(t) - 1,1);
        gamma(1) = 1/sqrt(1 - ((vz0^2 + vy0^2+vx0^2)/PhysConst.c^2));
        x=0.*t;
        x(1)=x0;
        y=0.*t;
        y(1)=y0;
        z=0.*t;
        z(1)=z0;
        pvz=0.*t;
        pvz(1)=gamma(1)*m*vz0;
        pvy=0.*t;
        pvy(1)=gamma(1)*m*vy0;
        pvx=0.*t;
        pvx(1)=gamma(1)*m*vx0;
        pvzd=0.*t;
        pvzd(1)=gamma(1)*m*vz0;
        pvyd=0.*t;
        pvyd(1)=gamma(1)*m*vy0;
        pvxd=0.*t;
        pvxd(1)=gamma(1)*m*vx0;
        vab=0.*t;
        vab(1) = vz0^2 + vy0^2 +vx0^2;
        vyd = 0.*t;
        vyd(1) = vy0;
        vzd = 0.*t;
        vzd(1) = vz0;
        vxd=0.*t;
        vxd(1)=vx0;
        gammad = zeros(numel(t) - 1,1);
        gammad(1) = 1/sqrt(1 - ((vz0^2 + vy0^2+vx0^2)/PhysConst.c^2));
        xd=0.*t;
        xd(1)=x0;
        yd=0.*t;
        yd(1)=y0;
        zd=0.*t;
        zd(1)=z0;
        vabd=0.*t;
        vabd(1) = vz0^2 + vy0^2 +vx0^2;
        eta=0.*t;
        eta(1)=omega*(t(1)-z(1)/c);
        envi=0.*t;
        
        envi(1) = -enveta(omega*(t(1)-z(1)/c),env);
        for i=2:numel(t)     
            eta(i) = omega*(t(i)-z(i-1)/c);
            envi(i) = -env.trans*enveta(eta(i),env);            
            pvz(i) =  pvz(1) - (0.5/(gamma(1)*m*(vz0-c)))*((e^2)*((envi(i)*sin(eta(i)))^2 -2*envi(i)*sin(eta(i))*envi(1)*sin(eta(1)) - (envi(1)*sin(eta(1)))^2) + 2*e*gamma(1)*m*vx0*envi(i)*sin(eta(i))); 
            pvx(i) = pvx(1) + e*(envi(i)*sin(eta(i)) - envi(1)*sin(eta(1)));
            pvy(i) = pvy(1);
            vab(i)=sqrt(((pvz(i)^2 + pvy(i)^2 + pvx(i)^2))/(c^2 * m^2 + (pvz(i)^2 + pvy(i)^2 + pvx(i)^2))); 
            gamma(i) = 1/(sqrt((1 + (vab(i)))*(1 - (vab(i)))));
            vx(i) = pvx(i)/(m*gamma(i));
            vz(i) = pvz(i)/(m*gamma(i));
            vy(i) = vy0;
            x(i) = x(i-1) + vx(i)*abs(t(i)-t(i-1));
            z(i) = z(i-1) + vz(i)*abs(t(i)-t(i-1));
            y(i) = y(i-1) + vy(i)*abs(t(i)-t(i-1));
            pvzd(i) = pvz(1) - (1/(2*gamma(1)*m*(vz(1)-c)))*(e^2 *(0.5*(envi(i)^2) + 2*envi(1)*sin(eta(1))*((envi(i)-envi(i-1))/(eta(i)-eta(i-1)))-(envi(1)*sin(eta(1)))^2) - 2*e*pvx(1)*((envi(i)-envi(i-1))/(eta(i)-eta(i-1))));
            pvyd(i) = pvy(1);
            pvxd(i) = pvx(1) - e*envi(1)*sin(eta(1)) + cos(eta(i))*e*((envi(i)-envi(i-1))/(eta(i)-eta(i-1)));
            vabd(i)=sqrt(((pvzd(i)^2 + pvyd(i)^2 + pvxd(i)^2))/(c^2 * m^2 + (pvzd(i)^2 + pvyd(i)^2 + pvxd(i)^2))); 
            gammad(i) = 1/(sqrt((1 + (vabd(i)))*(1 - (vabd(i)))));
            vxd(i) = pvxd(i)/(m*gammad(i));
            vzd(i) = pvzd(i)/(m*gammad(i));
            vyd(i) = vy0;
            xd(i) = xd(i-1) + vxd(i)*(t(i)-t(i-1));
            zd(i) = zd(i-1) + vzd(i)*(t(i)-t(i-1));
            yd(i) = yd(i-1) + vyd(i)*(t(i)-t(i-1));
            %disp((i/numel(t))*100)
        end    
    
end
input = struct('type', type, 'x0',x0,'y0',y0,'z0',z0, 'E0',E0, 'env',env,'vx0',vx0,'vy0',vy0,'vz0',vz0,'omega',omega,'t',t,'eta',eta);
info=struct('x',x,'y',y,'z',z,'vx',vx,'vy',vy,'vz',vz,'gamma',gamma,'xd',xd,'yd',yd,'zd',zd,'vxd',vxd,'vyd',vyd,'vzd',vzd,'gammad',gammad,'input',input,'pvx',pvx,'pvy',pvy,'pvz',pvz,'pvxd',pvxd,'pvyd',pvyd,'pvzd',pvzd,'vab',vab, 'ponder', pondermotive);
