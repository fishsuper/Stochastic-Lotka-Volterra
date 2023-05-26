function GradH = GradientH(y,a,b,gamma,mu)
syms y1 y2 y3 real
H = a*b*y1+y2-a*y3+gamma*log(y2)-mu*log(y3);
GradHy = jacobian(H,[y1,y2,y3]);
GradH = subs(GradHy,[y1,y2,y3],y);
