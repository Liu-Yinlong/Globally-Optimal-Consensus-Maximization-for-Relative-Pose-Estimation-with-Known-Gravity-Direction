function [X_3d] = exp_2d_2_3d(X_2d)
%EXP_2D_2_3D Summary of this function goes here
%   Detailed explanation goes here
rho=vecnorm(X_2d);
if(abs(rho)<1e-6)
    X_3d=[0;0;1];
    return ;
end
direc=X_2d./rho;

X_3d=[cos(rho);
    sin(rho).*direc];
end

