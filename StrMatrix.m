function B = StrMatrix(y,b,d)
B = [0,d*y(1)*y(2),b*d*y(1)*y(3);
    -d*y(1)*y(2),0,-y(2)*y(3);
    -b*d*y(1)*y(3),y(2)*y(3),0];