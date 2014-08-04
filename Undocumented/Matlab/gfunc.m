function avgint = gfunc(type, eta,envinf)
avgint=0.*eta;
etaac=0:2*pi/1000:2*pi;
if strcmp(type,'s')==1
    for i=1:numel(eta)
        for j=1:numel(etaac)
        avgint(i) = avgint(i) +  (1/(2*pi))*(enveta(eta(i)+etaac(j),envinf)*sin(etaac(j)))*(2*pi/1000);
        end
    end
elseif strcmp(type,'s2')==1
    for i=1:numel(eta)
        for j=1:numel(etaac)
        avgint(i) = avgint(i) +  (1/(2*pi))*(enveta(eta(i)+etaac(j),envinf)*(sin(etaac(j))^2))*(2*pi/1000);
        end
    end
elseif strcmp(type,'c')==1
    for i=1:numel(eta)
        for j=1:numel(etaac)
        avgint(i) = avgint(i) +  (1/(2*pi))*(enveta(eta(i)+etaac(j),envinf)*cos(etaac(j)))*(2*pi/1000);
        end
    end
elseif  stromp(type,'c2')==1
    for i=1:numel(eta)
        for j=1:numel(etaac)
        avgint(i) = avgint(i) +  (1/(2*pi))*(enveta(eta(i)+etaac(j),envinf)*(cos(etaac(j))^2))*(2*pi/1000);
        end
    end
end