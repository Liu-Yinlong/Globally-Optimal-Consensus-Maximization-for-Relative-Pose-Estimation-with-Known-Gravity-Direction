function [R_v,R_theta,t,x,y]=gen_data(num_inlier,num_outlier,noise_level)
%
%   t'*(corss(y',Rx))=0
%   R=exp([v*theta]x)
%
v=rand(3,1)*2-1;
R_v=v./norm(v);
if(R_v(1)<0)
    R_v=-R_v;
end
R_theta=rand()*2*pi-pi;
R=rotationVectorToMatrix(R_v*R_theta);

T=(rand(3,1)*2-1)*10;
t=T./norm(T);

X_3d=rand(3,num_inlier)*100+[500;500;500];

x_2d=X_3d./vecnorm(X_3d);
% T=[0;0;0];
Y_3d=R*X_3d+T+normrnd(0,noise_level,[3,num_inlier]);

y_2d=Y_3d./vecnorm(Y_3d);


x_outlier_non_norm=rand(3,num_outlier);
x_outlier=x_outlier_non_norm./vecnorm(x_outlier_non_norm);
x=[x_2d,x_outlier];

y_outlier_non_norm=rand(3,num_outlier);
y_outlier=y_outlier_non_norm./vecnorm(y_outlier_non_norm);
y=[y_2d,y_outlier];

end