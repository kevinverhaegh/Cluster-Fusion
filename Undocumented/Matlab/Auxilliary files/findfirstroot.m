function index=findfirstroot(x)
xsign=(abs(x)-x)./(2.*x);
for i=1:numel(x)
    if xsign(i+1) ~= xsign(i)
        index=i+1;
        return
    end
end
disp('No roots found!')
index=0;
return