clc;clear;close all

num_inlier=100;
num_outlier=10;

noise_level=0.0;

[R_v,R_theta,t_gt,x,y]=gen_data_5pt(num_inlier,num_outlier,noise_level);
t_=[0 -t_gt(3) t_gt(2);t_gt(3) 0 -t_gt(1);-t_gt(2) t_gt(1) 0];
R_gt=rotationVectorToMatrix(R_v*R_theta)
E_=t_*R_gt;
t_gt
[R_opt_out,t_opt_out,inlier_max] = ransac_5pt(x,y,R_v,0.0175*0.1);
% pts1=x./x(3,:);
% pts2=y./y(3,:);

% [E_opt, R_opt, t_opt, Eo] = five_point_algorithm(pts1(1:2,1:5), pts2(1:2,1:5), eye(3), eye(3));
% cameraParams = cameraParameters;
% [relativeOrientation,relativeLocation] = relativeCameraPose(E_,cameraParams,pts1(1:2,1:15)', pts2(1:2,1:15)');
% [rotationMatrix,translationVector] = cameraPoseToExtrinsics(relativeOrientation,relativeLocation);

