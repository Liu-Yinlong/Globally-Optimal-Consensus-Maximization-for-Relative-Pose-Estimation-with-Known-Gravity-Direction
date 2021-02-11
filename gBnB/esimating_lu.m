function [L,U,theta_opt] = esimating_lu(branch,a,b,c,epsilon)

t_c=exp_2d_2_3d(branch(1:2));
% t_r=branch(3);
[theta_opt,U]=get_opt_theta(branch,a,b,c,epsilon);%<- key part (inner BnB)

d=a+sin(theta_opt).*b+cos(theta_opt).*c;
% L_d=vecnorm(d);
% d_norm=d./repmat(L_d,3,1);

e=t_c'*d;
L=sum(abs(e)<=(epsilon));





end

