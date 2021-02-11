function [R_v,R_theta,t,x,y]=gen_data_bnb(num_inlier,num_outlier,noise_level)
%
%   t'*(corss(y',Rx))=0
%   R=exp([v*theta]x)
%
% v=rand(3,1)*2-1;
% R_v=v./norm(v);
% if(R_v(1)<0)
%     R_v=-R_v;
% end
% 
% 
% R_theta=(rand()*2-1)*pi/2;

R_v=[1;0;0];
R_theta=0;

R=rotationVectorToMatrix(R_v*R_theta);

num=num_inlier+num_outlier;

X_3d=(rand(3,num)*2-1)+[0;0;2];

X_3d_c=mean(X_3d,2);

T=-R*X_3d_c+X_3d_c+rand(3,1);
t=T./norm(T);

x=X_3d./vecnorm(X_3d);

noise_1=normrnd(0,noise_level,[2,num_inlier]);

noise_2=normrnd(0,20,[2,num_outlier]);

focal_length=1000;

Y_3d_=R*X_3d+T;
Y_3d_1=focal_length.*Y_3d_(1:2,1:num_inlier)./Y_3d_(3,1:num_inlier)+noise_1;

Y_3d_2=focal_length.*Y_3d_(1:2,end-num_outlier+1:end)./Y_3d_(3,end-num_outlier+1:end)+noise_2;

Y_3d=[[Y_3d_1,Y_3d_2];focal_length*ones(1,num)];
y=Y_3d./vecnorm(Y_3d);


end





