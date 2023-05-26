function F = EP_fun(y,x,stepsize,a,b,d,gamma,mu)

BGH = [(d*(x(1) + y(1))*(x(2) + y(2))*((gamma*(log((x(2))) - log((y(2)))))/(x(2) - y(2)) + 1))/4 - (b*d*(a + (mu*(log((x(3))) - log((y(3)))))/(x(3) - y(3)))*(x(1) + y(1))*(x(3) + y(3)))/4
    ((a + (mu*(log((x(3))) - log((y(3)))))/(x(3) - y(3)))*(x(2) + y(2))*(x(3) + y(3)))/4 - (a*b*d*(x(1) + y(1))*(x(2) + y(2)))/4
    ((x(2) + y(2))*(x(3) + y(3))*((gamma*(log((x(2))) - log((y(2)))))/(x(2) - y(2)) + 1))/4 - (a*b^2*d*(x(1) + y(1))*(x(3) + y(3)))/4];
F = y - x -BGH*stepsize;
