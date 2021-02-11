function [X_2d] = exp_3d_2_2d(X_3d)
%EXP_3D_2_2D Summary of this function goes here
%   Detailed explanation goes here
% Note: X_3d lies in upper hemi-sphere!
ind=X_3d(1,:)<0;
X_3d(:,ind)=-X_3d(:,ind);

rho=abs(acos(X_3d(1,:)));

dirc=X_3d(2:3,:)./sin(rho);
X_2d=rho.*dirc;

end

