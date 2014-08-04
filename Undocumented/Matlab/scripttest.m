for i=1:numel(avgint)
for j=1:numel(etaac)
%avgint(i) = avgint(i) +  (1/(2*pi))*exp(-(((eta(i)+etaac(j))^2)/(100^2)))*sin(etaac(j))*(2*pi/1000);
%avgint(i) = avgint(i) +  (1/(2*pi))*(cos((eta(i)+etaac(j))/(200*pi))^2)*sin(etaac(j))*(2*pi/1000);
avgint(i) = avgint(i) +  (1/(2*pi))*(((cos((eta(i)+etaac(j))/(200*pi))^2)*sin(etaac(j)))^2)*(2*pi/1000);
%avgint(i) = avgint(i) +  (1/(2*pi))*((exp(-(((eta(i)+etaac(j))^2)/(100^2)))*sin(etaac(j)))^2)*(2*pi/1000);
end
end