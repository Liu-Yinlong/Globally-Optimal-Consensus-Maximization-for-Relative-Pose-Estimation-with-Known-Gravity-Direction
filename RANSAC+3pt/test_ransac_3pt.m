clc;clear;close all

num_inlier=100;
num_outlier=1;
% num_outlier_list=[100:100:500];
% noise_level=0;
noise_level_list=0.001:0.001:0.01;
iter=100;
% outlier_num=length(num_outlier_list);
 noise_num=length(noise_level_list);
 for jj=1:noise_num
   

    for ii=1:iter
        [R_v,R_theta,t,x,y]=gen_data_3pt(num_inlier,num_outlier,noise_level_list(jj));

        x_unit=x./vecnorm(x);
        y_unit=y./vecnorm(y);
        epsilon=0.0175*0.5;
        [theta_opt,inlier_max(ii,jj),t_opt] = ransac_3pt(x_unit,y_unit,R_v,epsilon);

        angle_err(ii,jj)=abs(theta_opt-R_theta)*180/pi;
        t_err(ii,jj)=real(acos(abs(t'*t_opt)))*180/pi;
        disp([num2str(jj),'+',num2str(ii),'iter...'])
    end

end


