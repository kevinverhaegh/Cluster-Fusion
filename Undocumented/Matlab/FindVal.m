function yifra=FindVal(y,x,ifra)

xifra = (x(ceil(ifra)) - x(floor(ifra))) * (ifra - floor(ifra)) + x(floor(ifra));

a = (y(ceil(ifra)) - y(floor(ifra))) / (x(ceil(ifra)) - x(floor(ifra))) ;
b = y(ceil(ifra)) - a * x(ceil(ifra));

yifra = a * xifra + b;