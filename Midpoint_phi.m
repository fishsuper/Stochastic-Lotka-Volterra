function F = Mid_Phi(y,x,step)
global a b d gamma mu
BH = [d*(y(1)+x(1))/2*(gamma + (y(2)+x(2))/2) - b*d*(y(1)+x(1))/2*(a*(y(3)+x(3))/2 + mu)
      (y(2)+x(2))/2*((y(3)+x(3))/2*a + mu) - a*b*d*(y(1)+x(1))*(y(2)+x(2))/4
      - a*d*(y(1)+x(1))*(y(3)+x(3))*b^2/4 + (y(3)+x(3))/2*(gamma + (y(2)+x(2))/2)];
F = x + step*BH;