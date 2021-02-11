function [theta_opt,t_opt,inlier_max] = GBnB(x,y,R_v,epsilon)

% figure
% L_h=animatedline('color','r');
% U_h=animatedline('color','g');
% grid on
% legend('lower bound','upper bound')
% title('outer BnB convergence')
% xlabel('iteration')
% ylabel('inlier number')

%%%%%=========================
theta_opt=0;
inlier_max=0;

t_lu=[0;0;pi/2]; %center and radius

% x=x./sum(x.^2);
% y=y./sum(y.^2);

K=[0 -R_v(3) R_v(2);
    R_v(3) 0 -R_v(1);
    -R_v(2) R_v(1) 0]';

b=cross(y,K*x);
c=-cross(y,K^2*x);
a=cross(y,x)+cross(y,K^2*x);

%%%%%=========================

B=[];
L=0;
% U=size(x,2);

best_branch=t_lu;

iter=1;
t_opt_=[0;0];
while(1)
    sub=branch(best_branch);
    [L_1,U_1,theta_1]=esimating_lu(sub(:,1),a,b,c,epsilon);
    [L_2,U_2,theta_2]=esimating_lu(sub(:,2),a,b,c,epsilon);
    [L_3,U_3,theta_3]=esimating_lu(sub(:,3),a,b,c,epsilon);
    [L_4,U_4,theta_4]=esimating_lu(sub(:,4),a,b,c,epsilon);
       
    
    B=[B,[sub;[L_1 L_2,L_3,L_4;U_1,U_2,U_3,U_4;theta_1,theta_2,theta_3,theta_4]]];
    
    [L_temp,ind_L]=max(B(end-2,:));
    
    if(L_temp>L)
        L=L_temp;
        t_opt_=B(1:2,ind_L);
        theta_opt=B(end,ind_L);
        inlier_max=L;
    end
    
    temp_a=B(:,B(end-1,:)<L);
    
    B(:,B(end-1,:)<L)=[];
    
    
    [U,ind_U]=max(B(end-1,:));
    
    best_branch=B(1:3,ind_U);
    
    temp_b=B(:,ind_U);
    
    B(:,ind_U)=[];
    
%     addpoints(L_h,iter,L);
%     addpoints(U_h,iter,U);
%     drawnow
    
    if(L>=U)
        break;
    end
    
    iter=iter+1;
    
end

t_opt=exp_2d_2_3d(t_opt_);

end

%%%%%=========================

function new_b=branch(bra)

c=bra(1:2);
r=0.5*bra(3);

b1=[c(1)-r;c(2)-r;r];
b2=[c(1)+r;c(2)-r;r];
b3=[c(1)-r;c(2)+r;r];
b4=[c(1)+r;c(2)+r;r];

new_b=[b1 b2 b3 b4];

end