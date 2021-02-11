function [theta_opt,L] = get_opt_theta(branch,a,b,c,epsilon)
 
theta_opt=0;

t_c=exp_2d_2_3d(branch(1:2));
t_r=sqrt(2)*branch(3);

L=0;
% U=size(a,2);
Best_branch=[-pi;pi]/2;%<--------note this is the initial domain of R

iter=1;

B=[];

    while(1)
        center_b=0.5*(Best_branch(1)+Best_branch(2));
        new_b_1=[Best_branch(1);center_b ];
        new_b_2=[center_b; Best_branch(2)];

        [L_1,U_1]=get_LU_theta(new_b_1,a,b,c,epsilon,t_c,t_r);
        [L_2,U_2]=get_LU_theta(new_b_2,a,b,c,epsilon,t_c,t_r);

        B=[B [[new_b_1 new_b_2];[L_1,L_2;U_1 U_2]]];

        [L_temp,ind_L]=max(B(end-1,:));

        if(L_temp>L)
            L=L_temp;
            theta_opt=0.5*(B(1,ind_L)+B(2,ind_L));
        end

        B(:,B(end,:)<L)=[];

        [U,ind_U]=max(B(end,:));

        Best_branch=B(1:2,ind_U);

        B(:,ind_U)=[];

        iter=iter+1;
        if(L>=U)
            break;
        end

    end


end



function [L,U]=get_LU_theta(bran,a,b,c,epsilon,t_c,t_r)

center_b=0.5*(bran(1)+bran(2));
r_b=0.5*(bran(2)-bran(1));
d=a+sin(center_b).*b+cos(center_b).*c;
L_d=vecnorm(d);
d_norm=d./repmat(L_d,3,1);
% e=t_c'*d_norm;

item_1=sum(b.*b)*get_cos_2_max(bran(1),bran(2))+sum(c.*c)*get_sin_2_max(bran(1),bran(2));

[sin_cos_min,sin_cos_max]=get_sin_cos(bran(1),bran(2));

it_2=-2*dot(b,c);
f_2=it_2>0;
item_2=it_2*sin_cos_min;
item_2(f_2)=it_2(f_2)*sin_cos_max;

tau=sqrt(item_1+item_2)*r_b;


alpha_c=acos(t_c'*d_norm);
M=size(L_d,2);
alpha_min=max([zeros(1,M);alpha_c-t_r]);
alpha_max=min([pi*ones(1,M);alpha_c+t_r]);

L_min=L_d.*cos(alpha_max);
L_max=L_d.*cos(alpha_min);

final_min=L_min-tau;
final_max=L_max+tau;

flag_cross=final_min.*final_max<=0;
flag_pos=final_max(~flag_cross)>0 &final_min(~flag_cross)<=epsilon;
flag_neg=final_min(~flag_cross)<0 &abs(final_max(~flag_cross))<=epsilon;

U=sum(flag_cross)+sum(flag_pos)+sum(flag_neg);

flag_cross_lower=L_min.*L_max<=0;
flag_pos_lower=L_max(~flag_cross_lower)>0 &L_min(~flag_cross_lower)<=epsilon;
flag_neg_lower=L_min(~flag_cross_lower)<0 &abs(L_max(~flag_cross_lower))<=epsilon;

L=sum(flag_cross_lower)+sum(flag_pos_lower)+sum(flag_neg_lower);




end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=====================================
function [out_u]=get_sin_2_max(in_a,in_b)
flag_1=in_a<=-pi/2 &&-pi/2<=in_b;
flag_2=in_a<=pi/2 && pi/2<=in_b;
if(flag_1 || flag_2)
    out_u=1;
    return;
end

out_u=max(sin([in_a;in_b]).^2);

end
%=====================================
function [out_u]=get_cos_2_max(in_a,in_b)

flag_1=in_a<=0 &&0<=in_b;

if(flag_1)
    out_u=1;
    return;
end

out_u=max(cos([in_a;in_b]).^2);

end
%=====================================
function [out_l,out_u]=get_sin_cos(in_a,in_b)
a=2*in_a;
b=2*in_b;
flag_1=a<=pi/2 &&pi/2<=b;
flag_2=a<=-3*pi/2 &&-3*pi/2<=b;
flag_3=a<=-pi/2 &&-pi/2<=b;
flag_4=a<=3*pi/2 &&3*pi/2<=b;

out_u=0.5*max(sin([a;b]));
out_l=0.5*min(sin([a;b]));

if(flag_1||flag_2)
   out_u=0.5; 
end
if(flag_3||flag_4)
    out_l=-0.5;
end

end
