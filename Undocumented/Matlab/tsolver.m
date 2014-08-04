function Ztab=tsolver(N,t,R,F, alpha, Z)

Ztab = zeros(numel(t),numel(R));

for i=1:1:numel(t)
    for j=1:1:numel(R)
        [~, dum]=min(F(i)-((N*Z)/R(j)^2)*(1-3*(cos(alpha(:,j)./2)^2) + 2*cos(alpha(:,j)./2).^3));
        Ztab(i,j) = Z(dum);
    end
end
        