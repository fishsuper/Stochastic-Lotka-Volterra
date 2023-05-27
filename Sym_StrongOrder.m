%% setting paremeters of the equation %%
a = -2; b = -1; c = 0.5;
d = -0.5; gamma = 1; mu = 2;
y0 = [1.0;1.9;0.5];

%% setting attributes of the numerical scheme
T = 1; N = 2^12; 
dt = T/N; M = 500;
R0 = [1; 2; 4; 8; 16];

y_Mid = repmat(y0,1,M);   %% 3*M Matrix  
y_Sym = repmat(y0,1,M);
Y_Sym = zeros(3,M,5);
p_Sym = zeros(M,1); 
q_Sym = zeros(M,1); 
dW = sqrt(dt)*randn(M,N);  % M sample paths and N small intervals

stepsize = dt + c*dW;  % Midpoint rule
for j = 1:N
    for i = 1:M
        g = @(y)Mid_Fun(y,y_Mid(:,i),stepsize(i,j),a,b,d,gamma,mu);
        y_Mid(:,i) = fsolve(g,y_Mid(:,i));
    end
end

s = -log(y0(1))/d - b*log(y0(2)) + log(y0(3));
for p = 1:5
    p_Sym(:) = log(y0(1));
    q_Sym(:) = -log(y0(2))/d; % 每次循环都需初始化 
    R = R0(p); 
    Dt = R*dt; 
    L = N/R;
%     A = sqrt(Dt)*sqrt(4*abs(log(Dt)));
    for j = 1:L
        DetaW = sum(dW(:,(j-1)*R+1:j*R),2);
%         DetaW(DetaW>A) = A; DetaW(DetaW<-A) = -A;
        for i = 1:M
            g = @(p)Sym_P(a, b, d, s, p, q_Sym(i), gamma, mu, p_Sym(i), Dt, DetaW(i), c); % 关于p的函数，解p
            p_Sym(i) = fsolve(g,p_Sym(i));
            q_Sym(i) = Sym_Q(a, b, d, s, p_Sym(i), q_Sym(i), gamma, mu, Dt, DetaW(i), c);
        end
    end
    Y_Sym(1,:,p) = exp(p_Sym);
    Y_Sym(2,:,p) = exp(-d*q_Sym);
    Y_Sym(3,:,p) = exp(s+p_Sym/d-b*d*q_Sym);
end

D_SymNorm = zeros(M,5);
for p = 1:5
    for j = 1:M
        D_SymNorm(j,p) = norm(Y_Sym(:,j,p) - y_Mid(:,j));
    end
end
D_SymM = mean(D_SymNorm,1);

Dtvals = dt*R0;

h1 = loglog(Dtvals,D_SymM,'r*-'); hold on
h2 = loglog(Dtvals,Dtvals.^0.5,'b--');
xlabel('h','Fontsize',7)
title('MS errors','Fontsize',12)
set(gca,'Fontsize',7)
