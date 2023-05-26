function F = LP_phi(y,x,stepsize)
global I1 I2 I3

BGraH = [-((1/I2 - 1/I3)*(x(2) - y(2))*(x(3) - y(3)))/((log(x(2)) - log(y(2)))*(log(x(3)) - log(y(3))));
          ((1/I1 - 1/I3)*(x(1) - y(1))*(x(3) - y(3)))/((log(x(1)) - log(y(1)))*(log(x(3)) - log(y(3))));
         -((1/I1 - 1/I2)*(x(1) - y(1))*(x(2) - y(2)))/((log(x(1)) - log(y(1)))*(log(x(2)) - log(y(2))))];
F = x + BGraH*stepsize;