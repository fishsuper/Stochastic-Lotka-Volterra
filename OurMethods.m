% Wang Yuchao 2020-04-20
% compute stochastic L-V systems using our algorithm, Energy-preserving
% method and Midpoint scheme

global a b d gamma mu
a = -0.6; b = -1; c = 0;
d = -0.5; gamma = 1; mu = 2;
h = 0.001; T = 1000; t = 0:h:T;
y0 = [2.0;0.9;1.5];

N = length(t); M = 1;
y_IP = zeros(3,N,M); y_EP = zeros(3,N,M);
y_IP(:,1,:) = repmat(y0,1,M); y_EP(:,1,:) = repmat(y0,1,M);
y_Mid = zeros(3,N,M); y_Mid(:,1,:) = repmat(y0,1,M);
y_EM = zeros(3,N,M); y_EM(:,1,:) = repmat(y0,1,M);
normal = randn(N-1,M);
dW = sqrt(h)*normal;
Ah = sqrt(4*abs(log(h)));
normal(normal>Ah) = Ah; normal(normal<-Ah) = -Ah;
dW_hat = sqrt(h)*normal;
% Y1 = zeros(3,N-1,M);
% Z = zeros(N-1,M);
for j = 1:M
for i = 1:N-1
    stepsize = h+c*dW_hat(i,j);
    y1 = y_EP(:,i,j)+10; 
    y2 = y_EP(:,i,j)*(1+0.01*h); 
    y3 = y_Mid(:,i,j)+10; 
    y4 = y_Mid(:,i,j)+h;
    y5 = y_IP(:,i,j)+10;
    y6 = y_IP(:,i,j)*(1+0.01*h);
    while norm(y3-y4)>1e-8 && norm(y5-y6)>1e-8 && norm(y1-y2)>1e-8 
        y1 = y2;
        y2 = Energypreserving_phi(y1,y_EP(:,i,j),stepsize);
        y3 = y4;
        y4 = Midpoint_fun(y3,y_Mid(:,i,j),stepsize);
        y5 = y6;
        y6 = NewMethod_phi(y5,y_IP(:,i,j),stepsize);
%         y6 = Testf(y5,y_T(:,i,j),stepsize);
    end
    y_EP(:,i+1,j) = y2;
    y_Mid(:,i+1,j) = y4;
    y_IP(:,i+1,j) = y6;
    y_EM(:,i+1,j) = y_EM(:,i,j)+Euler_Increment(y_EM(:,i,j),h,dW(i,j),c);

end
end
H_Mid = Hamiltonian(y_Mid(:,:,1),a,b,gamma,mu);
H_EM = Hamiltonian(y_EM(:,:,1),a,b,gamma,mu);
H_EP = Hamiltonian(y_EP(:,:,1),a,b,gamma,mu);
H_IP = Hamiltonian(y_IP(:,:,1),a,b,gamma,mu);
C_Mid = -log(y_Mid(1,:,1))/d - b*log(y_Mid(2,:,1)) + log(y_Mid(3,:,1));
C_EM = -log(y_EM(1,:,1))/d - b*log(y_EM(2,:,1)) + log(y_EM(3,:,1));
C_EP = -log(y_EP(1,:,1))/d - b*log(y_EP(2,:,1)) + log(y_EP(3,:,1));
C_IP = -log(y_IP(1,:,1))/d - b*log(y_IP(2,:,1)) + log(y_IP(3,:,1));
k = 15000;
h1 = plot(t(1:k:end),C_IP(1:k:end),'k-o',t(1:k:end),C_EP(1:k:end),'r+-',t(1:k:end),C_Mid(1:k:end),'b-*',t(1:k:end),C_EM(1:k:end),'m-s');
set(h1,'LineWidth',1)
% ylim([-2,0.5])
legend('Our method','EP method','Midpoint method','Euler-Maruyama method')
title('Casimir','Fontsize',12)
ylim([-8,6])
ylabel('C(y_n)')
xlabel('t')

figure(2)
h2 = plot(t(1:k:end),H_IP(1:k:end),'k-o',t(1:k:end),H_EP(1:k:end),'r+-',t(1:k:end),H_Mid(1:k:end),'b-*',t(1:k:end),H_EM(1:k:end),'m-s');
set(h2,'LineWidth',1)
legend('Our method','EP method','Midpoint method','Euler-Maruyama method')
title('Hamiltonian','Fontsize',12)
ylabel('H(y_n)')
xlabel('t')
ylim([-6,4])


figure(3)
plot3(y_EM(1,:,1),y_EM(2,:,1),y_EM(3,:,1),'b')
xlabel('y^1')
ylabel('y^2')
zlabel('y^3')
title('Euler-Murayama method')
grid on

figure(4)
plot3(y_EP(1,:,1),y_EP(2,:,1),y_EP(3,:,1),'b')
xlabel('y^1')
ylabel('y^2')
zlabel('y^3')
title('EP method')
grid on

figure(5)
plot3(y_Mid(1,:,1),y_Mid(2,:,1),y_Mid(3,:,1),'b')
xlabel('y^1')
ylabel('y^2')
zlabel('y^3')
title('Midpoint method')
xlim([0,15])
zlim([0,20])
grid on

figure(6)
plot3(y_IP(1,:,1),y_IP(2,:,1),y_IP(3,:,1),'r')
xlabel('y^1')
ylabel('y^2')
zlabel('y^3')
title('Our method')
xlim([0,15])
zlim([0,20])
grid on

figure
subplot(131)
plot(t,y_IP(1,:,1),'Color','r','LineWidth',1)
hold on
plot(t,y_EP(1,:,1),'b--')
ylim([0,14])
xlim([0,20])
ylabel('y^1')
legend('Our method','EP method')
subplot(132)
plot(t,y_IP(2,:,1),'r','LineWidth',1)
hold on 
plot(t,y_EP(2,:,1),'b--')
ylim([0,1.6])
xlim([0,20])
ylabel('y^2')
subplot(133)
plot(t,y_IP(3,:,1),'r','LineWidth',1)
hold on
plot(t,y_EP(3,:,1),'b--')
xlim([0,20])
ylim([0,16])
ylabel('y^3')

%{
IP_M = mean(y_IP,3);
% EM_M = mean(y_EM,3);
H_IP = Hamiltonian(y_IP,a,b,gamma,mu);
H_EP = Hamiltonian(y_EP,a,b,gamma,mu);
C_IP = -log(y_IP(1,:,1))/d - b*log(y_IP(2,:,1)) + log(y_IP(3,:,1));
C_EP = -log(y_EP(1,:,1))/d - b*log(y_EP(2,:,1)) + log(y_EP(3,:,1));
% 
figure(4)
plot(t,C_IP,'k',t,C_EP,'b')
legend('Our proposed scheme','Energy preserving')
title('Carsimir')
figure(5)
plot(t,H_IP,'k',t,H_EP,'b')
legend('Our proposed scheme','Energy preserving')
title('Hamiltonian')
% 
figure(1)
hold on
subplot(322)
plot(t,y_Mid(1,:,1),'r',t,y_IP(1,:,1),'b')

% legend('Energy preserving scheme','Our proposed scheme')
ylabel('y1')
subplot(324)
plot(t,y_Mid(2,:,1),'r',t,y_IP(2,:,1),'b')

% legend('Energy preserving scheme','Our proposed scheme')
ylabel('y2')
subplot(326)
plot(t,y_Mid(3,:,1),'r',t,y_IP(3,:,1),'b')

% legend('Energy preserving scheme','Our proposed scheme')
ylabel('y3')

% subplot(322)
% h1 = plot(t,y_EP(1,:,1),'color',[.3 .3 .4]); hold on
% set(gca,'Fontsize',8)
% title('a = -0.2, b = -1, v = 0.5, \gamma = -1, \mu = 2')
% % set(h1,'Color','r','LineWidth',1);
% % plot(t,y_IP(1,:,1),'b--')
% ylabel('y_{1}')
% subplot(324)
% h1 = plot(t,y_EP(2,:,1),'color',[.3 .3 .4]); hold on
% set(gca,'Fontsize',8)
% ylabel('y_{2}')
% ylim([0,8])
% subplot(326)
% h1 = plot(t,y_EP(3,:,1),'color',[.3 .3 .4]); hold on
% set(gca,'Fontsize',8)
% ylim([0,50])
% ylabel('y_{3}')
% xlabel('(b)')

% subplot(321)
% h1 = plot(t,y_EP(1,:,1),'k'); hold on
% set(gca,'Fontsize',8)
% title('a = -2, b = -1, v = -0.5, \gamma = 1, \mu = 2')
% % set(h1,'Color','r','LineWidth',1);
% % plot(t,y_IP(1,:,1),'b--')
% ylabel('y_{1}')
% ylim([0,5])
% subplot(323)
% h1 = plot(t,y_EP(2,:,1),'k'); hold on
% set(gca,'Fontsize',8)
% ylabel('y_{2}')
% ylim([0,4])
% % ylim([0,8])
% subplot(325)
% h1 = plot(t,y_EP(3,:,1),'k'); hold on
% set(gca,'Fontsize',8)
% % ylim([0,50])
% ylabel('y_{3}')
% ylim([0,6])
% xlabel('(a)')
% Hamiltonian(y_EM,a,b,gamma,mu)
%}
%{
subplot(311)
plot(t,y_IP(1,:,1),'r--',t,y_EP(1,:,1),'b--',t,y_EM(1,:,1),'g')
title('y1')
legend('Invariant preserving','Energy preserving','Euler-Maruyama')
subplot(312)
plot(t,y_IP(2,:,1),'r--',t,y_EP(2,:,1),'b--',t,y_EM(2,:,1),'g')
title('y2')
subplot(313)
plot(t,y_IP(3,:,1),'r--',t,y_EP(3,:,1),'b--',t,y_EM(3,:,1),'g')
title('y3')
subplot(224)
plot(t,H_IP,'r',t,C_IP,'b')
title('Hamiltonian and Carsimir')
legend('Hamiltonian','Carsimir')
figure(2)
plot(t,C_IP,'r',t,C_EP,'b')
legend('Invariant preserving','Energy preserving')
title('Carsimir')
figure(3)
plot(t,H_IP,'r',t,H_EP,'b')
legend('Invariant preserving','Energy preserving')
title('Hamiltonian')
% Hamiltonian = a*b*y(1,:)+y(2,:)-a*y(3,:)+gamma*log(y(2,:))-mu*log(y(3,:))
%}
%}
