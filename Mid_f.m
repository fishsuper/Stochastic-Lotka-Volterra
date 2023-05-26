function F = Mid_f(y,x,stepsize)
BGH = [((x(1)+y(1))/2)*((x(2)+y(2))/2-2)
    ((x(2)+y(2))/2)*(-(x(1)+y(1))/2+1)];
F =x + BGH*stepsize;