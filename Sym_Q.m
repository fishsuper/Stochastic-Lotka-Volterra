function Q = Sym_Q(a, b, d, s, p, q, gamma, mu, h, zeta, c)
Kp = a*b*exp(p) - (a*exp(s + p/d - b*d*q))/d - mu/d;
KppKq = ((a*exp(s + p/d - b*d*q))/d^2 - a*b*exp(p))*(d*exp(-d*q) + d*(gamma - b*mu) - a*b*d*exp(s + p/d - b*d*q));
KpqKp = -a*b*exp(s + p/d - b*d*q)*(mu/d + (a*exp(s + p/d - b*d*q))/d - a*b*exp(p));
Q = q + 0.5*c^2*(KppKq+KpqKp)*h + Kp*(c*zeta+h);