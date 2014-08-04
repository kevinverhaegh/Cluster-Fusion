function [Z, IoDeg]=IonClusterChargePredict(input, Ef)
%detecting ion species

e0=1.602e-19;
mu0=4*pi*10^-7;
auE=5.14220652e11;
aut=2.418884326505e-17;
auEn=4.35974417e-18;
auL=5.2917720859e-11;
auv=2.1876912633e6;

IonSpecies=[0 0];
Ions=0;

for i=1:numel(input.AMUZcombo(:,1))
    if input.AMUZcombo(i,2)>0.01
        IonSpecies(end+1,:)=input.AMUZcombo(i,:);
        %counting end result ions in place
        Ions(end+1)=input.NumberInGroup(end, i);
    end
end

IonSpecies=IonSpecies(2:end, :);
Ions=Ions(2:end);

%Determine radius
R=max(input.x{2,1}(1,:));

QI=0;

%Determine Total Ion Charge
for i=1:numel(IonSpecies(:,1))
    QI=QI+IonSpecies(i,2)*Ions(i);
end

Z = (R/auL).*(R/auL).*(Ef/auE)/sqrt(2);

for i=1:numel(Z)
    if Z(i)>QI
        Z(i)=QI;
    end
end

for i=2:numel(Z)
    if Z(i-1)>Z(i)
        Z(i)=Z(i-1);
    end
end


IoDeg=Z./QI;