clc;clear;close all

addpath('./gBnB/')
addpath('./RANSAC+3pt')
addpath('./RANSAC+5pt')

num_all=100;

noise_level_list=0:0.5:2;

iter=50;
noise_num=length(noise_level_list);

outlier_rate=0.2;


epsilon=1e-3;
rho=0.99;

for jj=1:noise_num
     
        num_inlier=round(num_all*(1-outlier_rate));

        num_outlier=round(num_all-num_inlier);

        noise_level=noise_level_list(jj);
        
       ransac_3_iter_num=max(10,ceil(log(1-rho)./log(1-(1-outlier_rate).^3))+1);
       ransac_5_iter_num=max(10,ceil(log(1-rho)./log(1-(1-outlier_rate).^5))+1);

        
    for ii=1:iter
        
           
            [R_v,R_theta,t_gt,x,y]=gen_data(num_inlier,num_outlier,noise_level);
            R_gt=rotationVectorToMatrix(R_v*R_theta);
            x_unit=x./vecnorm(x);
            y_unit=y./vecnorm(y);
 
        %%%
        tic
        [theta_opt_3pt,t_opt_3pt,inlier_max_3pt(ii,jj)] = ransac_3pt(x_unit,y_unit,R_v,epsilon,ransac_3_iter_num);
        t_3pt(ii,jj)=toc;
        angle_err_3pt(ii,jj)=abs(((theta_opt_3pt-R_theta)))*180/pi;
        t_err_3pt(ii,jj)=real(acos(abs(t_gt'*t_opt_3pt)))*180/pi;
        
        %%%
        
        %%%
        tic
        [R_opt_5pt,t_opt_5pt,inlier_max_5pt(ii,jj)] = ransac_5pt(x_unit,y_unit,R_v,epsilon,ransac_5_iter_num);
        t_5pt(ii,jj)=toc;
        angle_err_5pt(ii,jj)=norm(rotationMatrixToVector(R_opt_5pt'*R_gt))*180/pi;
        t_err_5pt(ii,jj)=real(acos(abs(t_gt'*t_opt_5pt)))*180/pi;


        tic
        [theta_opt_bnb,t_opt_bnb,inlier_max_bnb(ii,jj)] = GBnB(x_unit,y_unit,R_v,epsilon);
        t_bnb(ii,jj)=toc;
        angle_err_bnb(ii,jj)=abs((theta_opt_bnb-R_theta))*180/pi;
        t_err_bnb(ii,jj)=real(acos(abs(t_gt'*t_opt_bnb)))*180/pi;
        
        disp([num2str(jj),'+',num2str(ii),'iter...'])
        
    end

end



%%

% figure
% hold on
% t_m_3pt=median(t_3pt);
% t_m_5pt=median(t_5pt);
% t_m_bnb=median(t_bnb);
% 
% plot(noise_level_list,t_m_3pt,'-r+')
% plot(noise_level_list,t_m_5pt,'-gs')
% plot(noise_level_list,t_m_bnb,'-b<')
% 
% legend(["RANSAC+3pt","RANSAC+5pt","gBnB"]);
% title("time")

figure
subplot(122)
hold on
t_err_a_3pt=mean(t_err_3pt);
t_err_a_5pt=mean(t_err_5pt);
t_err_a_bnb=mean(t_err_bnb);

plot(noise_level_list,t_err_a_3pt,'-r+')
plot(noise_level_list,t_err_a_5pt,'-gs')
plot(noise_level_list,t_err_a_bnb,'-k<')
legend(["RANSAC+3pt","RANSAC+5pt","gBnB"]);
title("Avg. trans error")
grid on
subplot(121)
hold on
r_err_a_3pt=mean(angle_err_3pt);
r_err_a_5pt=mean(angle_err_5pt);
r_err_a_bnb=mean(angle_err_bnb);
plot(noise_level_list,r_err_a_3pt,'-r+')
plot(noise_level_list,r_err_a_5pt,'-gs')
plot(noise_level_list,r_err_a_bnb,'-k<')
legend(["RANSAC+3pt","RANSAC+5pt","gBnB"]);
title("Avg. rot error")
grid on

figure
subplot(121)
hold on
in_m_3pt=median(inlier_max_3pt);
in_m_5pt=median(inlier_max_5pt);
in_m_bnb=median(inlier_max_bnb);

plot(noise_level_list,in_m_3pt,'-r+')
plot(noise_level_list,in_m_5pt,'-gs')
plot(noise_level_list,in_m_bnb,'-k<')

legend(["RANSAC+3pt","RANSAC+5pt","gBnB"]);
title("Median inlier number")
grid on

subplot(122)
hold on
success_inlier_threshold=2;

s_3pt=sum(t_err_3pt<=success_inlier_threshold & angle_err_3pt<=success_inlier_threshold)./iter;
s_5pt=sum(t_err_5pt<=success_inlier_threshold & angle_err_5pt<=success_inlier_threshold)./iter;
s_bnb=sum(t_err_bnb<=success_inlier_threshold & angle_err_bnb<=success_inlier_threshold)./iter;

plot(noise_level_list,s_3pt,'-r+')
plot(noise_level_list,s_5pt,'-gs')
plot(noise_level_list,s_bnb,'-k<')
legend(["RANSAC+3pt","RANSAC+5pt","gBnB"]);
title("success rate")
grid on















