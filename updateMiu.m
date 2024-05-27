function [miu,Z,e_miu] = updateMiu(F,lambda,gamma,Q,Z_old,miu_old,niu_old)
% 学习miu
[n,m] = size(F);
% Z = updateZ_miu_niu(F,lambda,gamma,Q,Z_old,miu_old,niu_old);
Z = Z_old;
clear Z_old

if m==n
    Y = 1 - sum(Z,2);
    e = n/100;
else
    Y = 1/n - sum(Z,2);
    e = 1/100;
end

iter = 0;
max_iter = 20;
while (iter <= max_iter) && ~(iter > 2 && sum(abs(Y)) < e)
    iter = iter + 1;
    d1Y = - sum( Z ./ (Z*gamma*2 + lambda) ,2);%一阶导数
    d2Y = sum( Z * lambda ./ (Z*gamma*2 + lambda).^3 ,2);%二阶导数
    
    % update miu
    miu = miu_old - (d1Y + sqrt(d1Y.^2 + Y.*d2Y*2)) ./ d2Y;% 二阶方法
    if ~isreal(sum(miu))
        miu = miu_old + Y ./ d1Y;% 一阶方法
    end
    miu_old = miu;
    Z = updateZ_miu_niu(F,lambda,gamma,Q,Z,miu_old,niu_old);
%     Z_old = Z;
    if m==n
        Y = 1 - sum(Z,2);
    else
        Y = 1/n - sum(Z,2);
    end
end

if m==n
    e_miu = sum(abs(sum(Z,2)-1)) + sum(abs(sum(Z,1)-1));
else
    e_miu = sum(abs(sum(Z,2)-1/n)) + sum(abs(sum(Z,1)-1/m));
end

end