function IncrementEM = Increment_EM(y,a,b,d,mu,gamma)
D1 = d*y(1).*(gamma + y(2)) - b*d*y(1).*(a*y(3) + mu);
D2 = y(2).*(y(3)*a + mu) - a*b*d*y(1).*y(2);
D3 = - a*d*y(1).*y(3)*b^2 +y(3).*(gamma + y(2));
IncrementEM = [D1;D2;D3];