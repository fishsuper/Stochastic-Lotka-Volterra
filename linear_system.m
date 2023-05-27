global B H1 H2
B = [0 1 -1; -1 0 3; 1 -3 0];
H1 = [2 1 1; 1 1 0; 1 0 1];
H2 = 0.25*[11 4 4;4 2 1;4 1 2];
h = 0.001; T = 200; t = 0:h:T;
y0 = [1 ; 0; -1];
N = length(t); M = 1;
y_mid = zeros(3,N,M); y_mid(:,1,:) = repmat(y0,1,M);
normal = randn(N-1,M);
dW = sqrt(h)*normal;
Ah = sqrt(4*abs(log(h)));
normal(normal>Ah) = Ah; normal(normal<-Ah) = -Ah;
dW_hat = sqrt(h)*normal;

for j = 1:M
for i = 1:N-1
    y1 = y_mid(:,i,j)+10; 
    y2 = y_mid(:,i,j)*(1+h); 
    while norm(y1-y2)>1e-8 
        y1 = y2;
        y2 = linear_midfun(y1,y_mid(:,i,j),h,dW_hat(i,j));
    end
    y_mid(:,i+1,j) = y2;
end
end
Cy_n = 3*y_mid(1,:,1) + y_mid(2,:,1) + y_mid(3,:,1);
plot(t, Cy_n)
Str = zeros(1,N);
I = eye(3);
for i = 1:N-1
    A = (I - 0.5*B*H1*h - 0.5*dW_hat(i,1).*B*H2)*B*(I - 0.5*B*H1*h - 0.5*dW_hat(i,1).*B*H2)'- (-I - 0.5*B*H1*h - 0.5*dW_hat(i,1).*B*H2)*B*(-I - 0.5*B*H1*h - 0.5*dW_hat(i,1).*B*H2)';
    Str(i) = norm(A);
end
plot(t,Str)
