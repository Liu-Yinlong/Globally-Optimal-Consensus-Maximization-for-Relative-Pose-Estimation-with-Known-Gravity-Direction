function [X,theta] = opt_3(K2,K1,K0)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
[X,E] = polyeig(K0,K1,K2);
ind=abs(imag(E))>1e-10;
E(ind)=[];
X(:,ind)=[];
theta=2*atan(real(E));
X=real(X);
end

