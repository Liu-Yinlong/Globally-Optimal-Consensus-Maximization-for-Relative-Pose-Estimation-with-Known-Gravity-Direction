clc;clear;close all

addpath('./gBnB_show/')

epsilon=1e-3;

num_inlier=100;
num_outlier=50;
noise_level=0.1;

[R_v,R_theta,t_gt,x,y]=gen_data(num_inlier,num_outlier,noise_level);
R_gt=rotationVectorToMatrix(R_v*R_theta);
x_unit=x./vecnorm(x);
y_unit=y./vecnorm(y);

tic
[theta_opt_bnb,t_opt_bnb,inlier_max_bnb] = GBnB_show(x_unit,y_unit,R_v,epsilon);
t_bnb=toc;

e_r=abs(theta_opt_bnb-R_theta)/pi*180;
e_t=acosd(abs(t_opt_bnb'*t_gt));

fprintf("time:%.2fs\t rot_err:%.2f(deg)\t tran_err:%.2f(deg)\n\n",t_bnb,e_r,e_t);





