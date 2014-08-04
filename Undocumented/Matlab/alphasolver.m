function alpha=alphasolver(N, Z, R, Ez)

alphaall=0:(1/100)*pi:pi;
alpha = zeros(numel(Z),numel(R));

for i=1:1:numel(Z)
    for j=1:1:numel(R)
        [~, dum]=min(abs(Ez(i).^2 ./ (4.*Z(i)) - ((N.*Z(i)) ./ (4.*R(j).^2)).*(2-3.*cos(alphaall) + cos(alphaall).^3 ).*(3 - sin(alphaall./2)).*(sin(alphaall./2).^2) + ((N.*Z(i))./R(j).^2).*(1-3.*cos(alphaall./2).^2 + 2.*(cos(alphaall./2).^3))));
        alpha(i,j) = alphaall(dum);
    end
end
        