function P = Sym_P(a, b, d, s, p, q, gamma, mu, pk, h, zeta, c)

Kq = a*b*d*exp(s + p/d - b*d*q) - d*(gamma - b*mu) - d*exp(-d*q);
KqpKq = -a*b*exp(s + p/d - b*d*q)*(d*exp(-d*q) + d*(gamma - b*mu) - a*b*d*exp(s + p/d - b*d*q));
KqqKp = -(d^2*exp(-d*q) - a*b^2*d^2*exp(s + p/d - b*d*q))*(mu/d + (a*exp(s + p/d - b*d*q))/d - a*b*exp(p));
P = p - pk + 0.5*c^2*(KqpKq+KqqKp)*h + (c*zeta+h)*Kq;
