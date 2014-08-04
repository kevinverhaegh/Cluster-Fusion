function [x,y,z]=IniGenSphere(R,N,rho,rhor)

try

rhoN=rho/max(rho);

x=zeros(N,1);
y=zeros(N,1);
z=zeros(N,1);

parfor i=1:N

    NotFound = 1;
    
    while NotFound==1
    xt=R*(2*rand-1);
    yt=R*(2*rand-1);
    zt=R*(2*rand-1);
    
    rt=sqrt(xt^2 + yt^2 + zt^2);
        
        if rt < R
        %Particle in sphere
        
        [~,ilow]=min(abs(rhor-rt));
        
        if ilow~= numel(rhor) && ilow ~= 1
        
        ihigh=ilow+1;
        
        deltar = rhor(ihigh)-rhor(ilow);
        
        ifra = (rt-rhor(ilow))/deltar + ilow;
        
        rhort = FindVal(rhoN,rhor,ifra);
        
        else
            
            if ilow == numel(rhor)
            
            ihigh = ilow;
            ilow = ihigh - 1;
            
            deltar = rhor(ihigh)-rhor(ilow);
        
            ifra = (rt-rhor(ilow))/deltar + ilow;
        
            rhort = FindVal(rhoN,rhor,ifra);
            
            else
                
                rhort = rhoN(1);
            end
        
        end
        
        if rand()<rhort
            %particle in distribution
            x(i)=xt;
            y(i)=yt;
            z(i)=zt;
            NotFound=0;
            
        end
            
        end
    end
end

catch
    disp(5)
end