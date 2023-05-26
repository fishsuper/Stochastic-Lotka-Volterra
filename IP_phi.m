function F = IP_phi(y,x,stepsize)
global a b d gamma mu
BGH = [(d*(x(1) - y(1))*(x(2) - y(2))*((gamma*(log(x(2)) - log(y(2))))/(x(2) - y(2)) + 1))/((log(x(1)) - log(y(1)))*(log(x(2)) - log(y(2)))) - (b*d*(a + (mu*(log(x(3)) - log(y(3))))/(x(3) - y(3)))*(x(1) - y(1))*(x(3) - y(3)))/((log(x(1)) - log(y(1)))*(log(x(3)) - log(y(3))))
    ((a + (mu*(log(x(3)) - log(y(3))))/(x(3) - y(3)))*(x(2) - y(2))*(x(3) - y(3)))/((log(x(2)) - log(y(2)))*(log(x(3)) - log(y(3)))) - (a*b*d*(x(1) - y(1))*(x(2) - y(2)))/((log(x(1)) - log(y(1)))*(log(x(2)) - log(y(2))))
    ((x(2) - y(2))*(x(3) - y(3))*((gamma*(log(x(2)) - log(y(2))))/(x(2) - y(2)) + 1))/((log(x(2)) - log(y(2)))*(log(x(3)) - log(y(3)))) - (a*b^2*d*(x(1) - y(1))*(x(3) - y(3)))/((log(x(1)) - log(y(1)))*(log(x(3)) - log(y(3))))];
F = x + BGH*stepsize;