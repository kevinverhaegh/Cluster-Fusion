function [specfun, iflag] = f_solve_obj_NAG(n,x,iflag)
global Qi s0i taui
specfun=zeros(n,1);
for k=1:n
specfun(k)=( (sqrt(x(k)*(x(k)-1)) + log(sqrt(x(k)) + sqrt(x(k)-1)) - taui*sqrt(2*Qi/(3*s0i^3)))) + (x(k)>1);
end
aap=1;