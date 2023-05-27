function F = midfunction(y,x,h,Wh)
global B H1 H2
F =x + B*(H1*(x+y)*h/2 + H2*(x+y)*Wh/2);
