function [theta_opt,t_opt,inlier_max] = ransac_3pt(x,y,R_v,epsilon,ransac_3_iter_num)
%RANSAC_3PT 此处显示有关此函数的摘要
%   此处显示详细说明

x=x./vecnorm(x);
y=y./vecnorm(y);

K=[0 -R_v(3) R_v(2);
    R_v(3) 0 -R_v(1);
    -R_v(2) R_v(1) 0]';

B=cross(y,K*x);
C=-cross(y,K^2*x);
A=cross(y,x)+cross(y,K^2*x);

K0=A+C;
K1=2*B;
K2=A-C;

num_iter=ransac_3_iter_num;
N=size(K0,2);
ind=randi(N,3,num_iter);

theta_list=[];
inlier_num_list=[];
t_list=[];

t_start=tic;


for ii=1:num_iter
    ind_=ind(:,ii);
    [t_,theta_opt_] = opt_3(K2(:,ind_')',K1(:,ind_')',K0(:,ind_')');
    num_theta=length(theta_opt_);
    inlier_num=zeros(num_theta,1);
    for jj=1:num_theta
        d=(A+sin(theta_opt_(jj)).*B+cos(theta_opt_(jj)).*C);
        angle_err=t_(:,jj)'*d;
        inlier_num(jj)=sum(abs((angle_err))<=(epsilon));
    end
    inlier_num_list=[inlier_num_list;inlier_num];
    theta_list=[theta_list;theta_opt_];
    t_list=[t_list,t_];
    
    
    t_end=toc(t_start);
    
    if(t_end>=60)%% 1min=60 seconds
        break;
    end
    
end

[inlier_max,ind_max]=max(inlier_num_list);
theta_opt=theta_list(ind_max);
t_opt=t_list(:,ind_max);
end

