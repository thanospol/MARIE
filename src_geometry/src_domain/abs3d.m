function absR = abs3d(r)
absR2 = abs(r).^2;
absR = sqrt(sum(absR2, 4));