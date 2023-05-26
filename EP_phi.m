function F = EP_phi(y,x,stepsize)
% x = y0;
% y = y0+h;
global a b d gamma mu
BGH = [(d*(x(1) + y(1))*(x(2) + y(2))*((gamma*(log(abs((x(2)))) - log(abs((y(2))))))/(x(2) - y(2)) + 1))/4 - (b*d*(a + (mu*(log(abs((x(3)))) - log(abs((y(3))))))/(x(3) - y(3)))*(x(1) + y(1))*(x(3) + y(3)))/4
    ((a + (mu*(log(abs((x(3)))) - log(abs((y(3))))))/(x(3) - y(3)))*(x(2) + y(2))*(x(3) + y(3)))/4 - (a*b*d*(x(1) + y(1))*(x(2) + y(2)))/4
    ((x(2) + y(2))*(x(3) + y(3))*((gamma*(log(abs((x(2)))) - log(abs((y(2))))))/(x(2) - y(2)) + 1))/4 - (a*b^2*d*(x(1) + y(1))*(x(3) + y(3)))/4];
F = x + BGH*stepsize;