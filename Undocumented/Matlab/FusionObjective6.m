function optim = FusionObjective6(x)

input.R=x(1)*1e-7;
input.n0=x(2)*1e26;
input.RhoA = x(3);
input.RhoB = x(4);

output=main_file(input);

optim=output.Ej / (2.45 * 1e6 * PhysConst.e * output.Nfusion);

disp('Q-factor:')
disp(1/optim);

dlmwrite('FusionObjective6Results',horzcat(input.R,input.n0,input.RhoA,input.RhoB,output.Nfusion,output.expansion.FusY,output.Ej,output.npack),'delimiter','\t','-append')


