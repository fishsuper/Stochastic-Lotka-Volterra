function H = Hamiltonian(y,a,b,gamma,mu)
H = a*b*y(1,:)+y(2,:)-a*y(3,:)+gamma*log(y(2,:))-mu*log(y(3,:));