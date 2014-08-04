function optim = FusionObjective5(x)

input.R=x(1);
input.n0=x(2);
input.RhoA = x(3);
input.RhoB = x(4);

output=main_file(input);

optim=1 / output.Nfusion;

dlmwrite('FusionObjective5Results',horzcat(input.R,input.n0,input.RhoA,input.RhoB,output.Nfusion,output.expansion.FusY,output.Ej,output.npack),'delimiter','\t','-append')


