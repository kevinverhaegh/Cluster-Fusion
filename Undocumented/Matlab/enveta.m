function env=enveta(eta,envinf)

if strcmp('gaussian',envinf.type)==1
    env = (envinf.A/envinf.omega) .* exp(-(eta./envinf.omega).^2 ./ (2 .* envinf.tau .* envinf.tau));
    
elseif strcmp('cos2',envinf.type)==1
    env = envinf.A * (cos((eta/(omega*(envinf.dura/(2*pi)))))^2);
    
elseif strcmp('tri',envinf.type)==1
    
elseif strcmp('square',envinf.type)==1
    
end