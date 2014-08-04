function Pres = Bayes(P0, Pmat)

Pres = zeros(numel(Pmat(:,1))+1,numel(Pmat(1,:)));
Pres(1,:) = P0;
for i=1:numel(Pmat(:,1))
    Pres(i+1,1)=(Pmat(i,1)*Pres(i,1))/(Pmat(i,1)*Pres(i,1) + Pmat(i,2)*Pres(i,2));
    Pres(i+1,2)=(Pmat(i,2)*Pres(i,2))/(Pmat(i,1)*Pres(i,1) + Pmat(i,2)*Pres(i,2));
end