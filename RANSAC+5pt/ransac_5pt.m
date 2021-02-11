function [R_opt_out,t_opt_out,inlier_max] = ransac_5pt(x,y,R_v,epsilon,ransac_5_iter_num)
%RANSAC_5PT 此处显示有关此函数的摘要
%   此处显示详细说明
R_list=[];
inlier_num_list=[];
t_list=[];

num_iter=ransac_5_iter_num;
N=size(x,2);
ind=randi(N,5,num_iter);

t_start=tic;

for ii=1:num_iter
    ind_=ind(:,ii);
    
    pts1=x(1:2,ind_)./x(3,ind_);
    pts2=y(1:2,ind_)./y(3,ind_);

    [E_opt, R_opt, t_opt, ~] = five_point_algorithm(pts1, pts2, eye(3), eye(3));
    M=size(E_opt,1);
    inlier_num=zeros(M,1);
    for jj=1:M
        d=cross(y,R_opt{jj}*x);
        angle_err=t_opt{jj}'*d;
        inlier_num(jj)=sum(abs((angle_err))<=(epsilon));
    end
    
    inlier_num_list=[inlier_num_list;inlier_num];
    R_list=[R_list;R_opt];
    t_list=[t_list;t_opt];
    
    t_end=toc(t_start);
    
    if(t_end>=60)%% 1min=60 seconds
        break;
    end
    
end

    [inlier_max,ind_max]=max(inlier_num_list);
    R_opt_out=R_list{ind_max};
    t_opt_out=t_list{ind_max};

end

