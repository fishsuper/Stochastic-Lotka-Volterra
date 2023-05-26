function F = EP_f(y,x,stepsize)
BGH = [-0.25*(x(1) + y(1))*(x(2) + y(2))*(-1 - 0.5*(log(abs(y(2))) - log(abs(x(2))))/(y(2) - x(2)))
0.25*(x(1) + y(1))*(x(2) + y(2))*(0.4 + (log(abs(y(1))) - log(abs(x(1))))/(y(1) - x(1)))];
F =x + BGH*stepsize;