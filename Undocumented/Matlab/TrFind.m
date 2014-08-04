function F = TrFind(x, pred1, pred2, t, i)
tr = t(i) - sqrt(x(1)^2 + x(2)^2 + x(3)^2)/PhysConst.c;
F = [pred1.x(i) - evalu(pred2.x, t, tr) - x(1);
    pred1.y(i) - evalu(pred2.y, t, tr) - x(2);
    pred1.z(i) - evalu(pred2.z, t, tr) - x(3)];
    