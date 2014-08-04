function Te = TeDistr(Bx,By,Bz)


    Lorentz = 1./sqrt(1-(Bx.^2 + By.^2 + Bz.^2));
    Te = (((2/3).* PhysConst.me .* PhysConst.c^2 .* (Lorentz - 1))./(PhysConst.e));
