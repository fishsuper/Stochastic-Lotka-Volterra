%% setting paremeters of the equation %%
a = -0.6; b = -1; c = 0;
d = -0.5; gamma = 1; mu = 2;
y0 = [2.0;0.9;1.5];

%% setting attributes of the numerical scheme
T = 1; N = 2^12; 
dt = T/N; M = 1000;
R0 = 2*[1; 2; 4; 8; 16];
y_EP = repmat(y0,1,M);
Y_EP = zeros(3,M,5);
y_EM = repmat(y0,1,M);
Y_EM = zeros(3,M,5);

y_Mid = repmat(y0,1,M);   %% 3*M Matrix  
y_IP = repmat(y0,1,M);
Y_IP = zeros(3,M,5);
% M sample paths and N small intervals
normal = randn(M,N);
W = sqrt(dt)*normal;
Ah = sqrt(4*abs(log(dt)));
normal(normal>Ah) = Ah; normal(normal<-Ah) = -Ah;
dW = sqrt(dt)*normal;

stepsize = dt + c*dW;  % Midpoint rule
for j = 1:N
    for i = 1:M
        y = y_Mid(:,i)+10; 
        y2 = y_Mid(:,i)+dt; 
        x = y_Mid(:,i);
        while norm(y-y2) > 1e-8
            y = y2;
            BH = [d*(y(1)+x(1))/2*(gamma + (y(2)+x(2))/2) - b*d*(y(1)+x(1))/2*(a*(y(3)+x(3))/2 + mu)
                   (y(2)+x(2))/2*((y(3)+x(3))/2*a + mu) - a*b*d*(y(1)+x(1))*(y(2)+x(2))/4
                    - a*d*(y(1)+x(1))*(y(3)+x(3))*b^2/4 + (y(3)+x(3))/2*(gamma + (y(2)+x(2))/2)];
            y2 = x + stepsize(i,j)*BH;
%             y = y2;
%             y2 = Mid_Phi(y,y_Mid(:,i),stepsize(i,j));
        end
        y_Mid(:,i) = y2;
    end
end

for p = 1:5
    y_IP = repmat(y0,1,M);
    y_EP = repmat(y0,1,M);
    R = R0(p); 
    Dt = R*dt; 
%     sqrtDt = sqrt(Dt);
    L = N/R;
%     A = sqrt(Dt)*sqrt(4*abs(log(Dt)));
    for j = 1:L
        DetaW = sum(dW(:,(j-1)*R+1:j*R),2);
%         DetaW(DetaW>A) = A; DetaW(DetaW<-A) = -A;
        for i = 1:M
            y = y_IP(:,i)+10; 
            y4 = y_IP(:,i)*(1+0.01*Dt); 
            yy = y_EP(:,i)+10; 
            y6 = y_EP(:,i)*(1+0.01*Dt);
            x = y_IP(:,i);
            xx = y_EP(:,i);
            while norm(y-y4)>1e-8 && norm(yy-y6)>1e-8
                y = y4;
                BGH = [(d*(x(1) - y(1))*(x(2) - y(2))*((gamma*(log(abs(x(2))) - log(abs(y(2)))))/(x(2) - y(2)) + 1))/((log(abs(x(1))) - log(abs(y(1))))*(log(abs(x(2))) - log(abs(y(2))))) - (b*d*(a + (mu*(log(abs(x(3))) - log(abs(y(3)))))/(x(3) - y(3)))*(x(1) - y(1))*(x(3) - y(3)))/((log(abs(x(1))) - log(abs(y(1))))*(log(abs(x(3))) - log(abs(y(3)))))
                       ((a + (mu*(log(abs(x(3))) - log(abs(y(3)))))/(x(3) - y(3)))*(x(2) - y(2))*(x(3) - y(3)))/((log(abs(x(2))) - log(abs(y(2))))*(log(abs(x(3))) - log(abs(y(3))))) - (a*b*d*(x(1) - y(1))*(x(2) - y(2)))/((log(abs(x(1))) - log(abs(y(1))))*(log(abs(x(2))) - log(abs(y(2)))))
                        ((x(2) - y(2))*(x(3) - y(3))*((gamma*(log(abs(x(2))) - log(abs(y(2)))))/(x(2) - y(2)) + 1))/((log(abs(x(2))) - log(abs(y(2))))*(log(abs(x(3))) - log(abs(y(3))))) - (a*b^2*d*(x(1) - y(1))*(x(3) - y(3)))/((log(x(1)) - log(abs(y(1))))*(log(abs(x(3))) - log(abs(y(3)))))];
                y4 = x + BGH*(Dt+c*DetaW(i));
                yy = y6;
                EBGH = [(d*(xx(1) + yy(1))*(xx(2) + yy(2))*((gamma*(log(abs((xx(2)))) - log(abs((yy(2))))))/(xx(2) - yy(2)) + 1))/4 - (b*d*(a + (mu*(log(abs((xx(3)))) - log(abs((yy(3))))))/(xx(3) - yy(3)))*(xx(1) + yy(1))*(xx(3) + yy(3)))/4
                    ((a + (mu*(log(abs((xx(3)))) - log(abs((yy(3))))))/(xx(3) - yy(3)))*(xx(2) + yy(2))*(xx(3) + yy(3)))/4 - (a*b*d*(xx(1) + yy(1))*(xx(2) + yy(2)))/4
                    ((xx(2) + yy(2))*(xx(3) + yy(3))*((gamma*(log(abs((xx(2)))) - log(abs((yy(2))))))/(xx(2) - yy(2)) + 1))/4 - (a*b^2*d*(xx(1) + yy(1))*(xx(3) + yy(3)))/4];
                y6 = xx + EBGH*(Dt+c*DetaW(i));
            end
            y_IP(:,i) = y4;
        end
    end
    Y_IP(:,:,p) = y_IP;
    Y_EP(:,:,p) = y_EP;
    Y_EM(:,:,p) = y_EM;
end

D_IPNorm = zeros(M,5);
D_EPNorm = zeros(M,5);
D_EMNorm = zeros(M,5);
for p = 1:5
    for j = 1:M
        D_IPNorm(j,p) = norm(Y_IP(:,j,p) - y_Mid(:,j));
        D_EPNorm(j,p) = norm(Y_EP(:,j,p) - y_Mid(:,j));
        D_EMNorm(j,p) = norm(Y_EM(:,j,p) - y_Mid(:,j));
    end
end
D_IPM = mean(D_IPNorm,1);
D_EPM = mean(D_EPNorm,1);
D_EMM = mean(D_EMNorm,1);

Dtvals = dt*R0;
figure
h1 = loglog(Dtvals,D_IPM,'ko-'); hold on
h2 = loglog(Dtvals,D_EPM,'r+-'); hold on
set(h1,'linewidth',1)
set(h2,'linewidth',1)
% h4 = loglog(Dtvals,D_EMM,'r+-'); hold on
h3 = loglog(Dtvals,Dtvals.^2,'b--');
% loglog(Dtvals,Dtvals.^0.5,'r--');
xlabel('h','Fontsize',8)
title('MS errors','Fontsize',12)
set(gca,'Fontsize',8)
legend('Our method','EP method','line of slope 2')