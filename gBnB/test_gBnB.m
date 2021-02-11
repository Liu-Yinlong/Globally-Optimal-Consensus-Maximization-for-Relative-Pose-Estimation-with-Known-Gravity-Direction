clc;clear;close all

num_inlier=50;
num_outlier=20;
noise_level=0.5;
for ii=1:10
    
[R_v,R_theta,t_gt,x,y]=gen_data_bnb(num_inlier,num_outlier,noise_level);
epsilon=1e-3;
tic
[theta_opt,t_opt,inlier_num(ii)] = GBnB(x,y,R_v,epsilon);
tim(ii)=toc;
e_r(ii)=abs(theta_opt-R_theta)/pi*180;
e_t(ii)=acos(abs(t_opt'*t_gt))/pi*180;


R_gt=rotationVectorToMatrix(R_v*R_theta);
error_E=abs(t_gt'*cross(y,R_gt*x));
Inlier_gt(ii)=sum((error_E)<=epsilon);

R_opt=rotationVectorToMatrix(R_v*theta_opt);
error_opt=abs(t_opt'*cross(y,R_opt*x));

disp([num2str(ii),' iter...']);

if(e_t(ii)>5)
    keyboard
end

end

disp(['R error(deg):' num2str(e_r)])
disp(['T error(deg):' num2str(e_t)])
figure
subplot(121)
bar(e_r)
ylabel('rot err')
subplot(122)
bar(e_t)
ylabel('tran err')
figure
plot(inlier_num,'ro-');
hold on
plot(Inlier_gt,'gs-');
ylabel('inlier num');
title('global solution Vs. ground truth');
legend('global solution','ground truth');
