function [Z,miu,niu] = updateZ(F,lambda,gamma,Q,Z_old,miu_old,niu_old)
% min F.*Z + lambda*Z.*logZ + gamma*Z.^2
[n,m] = size(F);
if nargin < 4
    if m==n
        Q = ones(n,n) - eye(n);%方阵对角线为0
    else
        Q = ones(n,m);%全1
    end
end
Q=(Q==1);%转为逻辑值logical，计算量为一次运算
if nargin < 5
    if m==n
        Z_old = Q./repmat(sum(Q,2),1,m);%方阵行和为1
    else
        Z_old = Q./repmat(sum(Q,2),1,m)/n;%行和为1/n
    end
end
if nargin < 6 || isempty(miu_old)
    Temp = F + (log(Z_old)+1)*lambda + Z_old*gamma*2;
    Temp(~Q) = 0;
%     sum_F = - sum(Temp(Q)) /(m*n*2);
%     miu_old = sum_F *ones(n,1);%用Z平均初始化乘子
%     niu_old = sum_F *ones(m,1);
    miu_old = - sum(Temp,2) /(m*2);%用Z平均初始化乘子/3
    niu_old = - sum(Temp)' /(n*2);
    clear Temp
%     Z_temp = zeros(n,m);
%     Temp = F + (log(Z_old)+1)*lambda + Z_old*gamma*2;
%     Z_temp(Q) = Temp(Q);
%     sum_F = - sum(Z_temp,'all') /(m*n*2);
%     miu_old = sum_F *ones(n,1);%用Z平均初始化乘子
%     niu_old = sum_F *ones(m,1);
%     clear Z_temp Temp
end

Z = updateZ_miu_niu(F,lambda,gamma,Q,Z_old,miu_old,niu_old);
miu = miu_old;
niu = niu_old;
% Z_old = Z;
clear Z_old miu_old niu_old
if m==n
    e = n/100;
else
    e = 1/100;
end

iter = 0;
max_iter = 20;
while (iter <= max_iter) && ~(iter > 2 && e_miu + e_niu < e)
    iter = iter + 1;
    [miu,Z,e_miu] = updateMiu(F,lambda,gamma,Q,Z,miu,niu);
    [niu,Z,e_niu] = updateNiu(F,lambda,gamma,Q,Z,miu,niu);
end
%     [e_miu + e_niu e]
% 1

