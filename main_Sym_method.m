%%%%%% Wang YC 2020-04-20
%%%%%% solving the stochastic Lotka-Volterra systems via symplectic method

%% inatial paremeter-values of the equation
a = -2; b = -1; c = 0.5;
d = -0.5; gamma = 1; mu = 2;
y0 = [1.0;1.9;0.5];
%%
s = -log(y0(1))/d - b*log(y0(2)) + log(y0(3));
h = 0.001; Ah = sqrt(4*abs(log(h)));
T = 0.5; t = 0:h:T;
n = length(t); m = 1;
Normal = randn(n-1,m);
% dW = sqrt(h)*Normal;
Normal(Normal>Ah) = Ah; Normal(Normal<-Ah) = -Ah;
dW_hat = sqrt(h)*Normal;
y_Sym = zeros(3,n,m);
y_Sym(:,1,:) = repmat(y0,1,m);
p_Sym = zeros(n,m); p_Sym(1,:) = log(y0(1))/b/d;
q_Sym = zeros(n,m); q_Sym(1,:) = log(y0(3));
%s pk
%% symplectic method
for j = 1:m
    for i = 1:n-1
        g = @(p)Sym_P(a, b, d, s, p, q_Sym(i,j), gamma, mu, p_Sym(i,j), h, dW_hat(i,j), c);
        p_Sym(i+1,j) = fsolve(g,p_Sym(i,j));
        q_Sym(i+1,j) = Sym_Q(a, b, d, s, p_Sym(i+1,j), q_Sym(i,j), gamma, mu, h, dW_hat(i,j), c);
    end
end
y_Sym(1,:,:) = exp(b*d*p_Sym);
y_Sym(3,:,:) = exp(q_Sym);
y_Sym(2,:,:) = exp((q_Sym-s)/b-p_Sym);
Sym_M = mean(y_Sym,3);
for k = 1:3
    figure(k)
    plot(t,Sym_M(k,:),'k')
end
