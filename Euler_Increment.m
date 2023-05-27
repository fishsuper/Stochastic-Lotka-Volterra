function inc = Euler_Increment(y,Dt,dW,c)
global a b d gamma mu
% y(1,:) = y(1);
% y(2,:) = y(2);
% y(3,:) = y(3);
BGH = [ y(1,:)*(d*gamma + d*y(2,:) - b*d*mu - a*b*d*y(3,:))
           y(2,:)*(mu + a*y(3,:) - a*b*d*y(1,:))
         y(3,:)*(- a*d*y(1,:)*b^2 + gamma + y(2,:))];
PBGH = [y(1,:)*(d*gamma + d*y(2,:) - b*d*mu - a*b*d*y(3,:))^2 + d*y(1,:)*y(2,:)*(mu + a*y(3,:) - a*b*d*y(1,:)) - a*b*d*y(1,:)*y(3,:)*(- a*d*y(1,:)*b^2 + gamma + y(2,:))
        y(2,:)*(mu + a*y(3,:) - a*b*d*y(1,:))^2 + a*y(2,:)*y(3,:)*(- a*d*y(1,:)*b^2 + gamma + y(2,:)) - a*b*d*y(1,:)*y(2,:)*(d*gamma + d*y(2,:) - b*d*mu - a*b*d*y(3,:))
        y(3,:)*(- a*d*y(1,:)*b^2 + gamma + y(2,:))^2 + y(2,:)*y(3,:)*(mu + a*y(3,:) - a*b*d*y(1,:)) - a*b^2*d*y(1,:)*y(3,:)*(d*gamma + d*y(2,:) - b*d*mu - a*b*d*y(3,:))];
inc = BGH*(Dt+c*dW)+Dt*c^2/2*PBGH;