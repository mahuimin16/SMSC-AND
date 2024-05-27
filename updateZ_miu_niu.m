function Z = updateZ_miu_niu(F,lambda,gamma,Q,Z_old,miu,niu)
% 固定miu_niu学习Z
[n,m] = size(F);
if m==n
    e = 1/n/100;%定义一个小量
else
    e = 1/n/m/100;
end
% Z = Z_old;
% clear Z_old

% for i=1:30 % 一阶方法
%     Z = -Z.*( F + lambda*log(Z) + repmat(miu,1,m) + repmat(niu',n,1)) ./ (Z.*gamma*2 + lambda);
%     Z = (abs(Z - e) + Z - e) / 2 + e;%让所有Z都>e，Q=0对应的是nan
% %     Z_old = Z;
% end

% for i=1:3 % 二阶方法
%     Y = F + lambda*(log(Z)+1) + Z.*gamma*2 + repmat(miu,1,m) + repmat(niu',n,1);
%     Z = Z*2 + Z.*( Z.*gamma*2 - sqrt((Z.*gamma*2+lambda).^2 + Y*lambda*2) ) ./ lambda;
% %     Z_old = Z;
% end

for i=1:10% 一阶方法，避免陷入循环
    Z = real(-Z_old.*( F + lambda*log(Z_old) + repmat(miu,1,m) + repmat(niu',n,1)) ./ (Z_old.*gamma*2 + lambda));
    Z = (abs(Z - e) + Z - e) / 2 + e;%让所有Z都>e，Q=0对应的是nan
    Z_old = Z;
end

for i=1:15 % 二阶方法
    Y = real(F + lambda*(log(Z_old)+1) + Z_old.*gamma*2 + repmat(miu,1,m) + repmat(niu',n,1));
    Z = Z_old*2 + Z_old.*( Z_old.*gamma*2 - sqrt((Z_old.*gamma*2+lambda).^2 + Y*lambda*2) ) ./ lambda;
    %     Z_old = Z;
%     ~isreal(Z)
    if ~isreal(sum(Z,'all')) % 有复数，用一阶方法
        Z = real(-Z_old.*( F + lambda*log(Z_old) + repmat(miu,1,m) + repmat(niu',n,1)) ./ (Z_old.*gamma*2 + lambda));
%         ~isreal(Z)
    end
    Z = (abs(Z - e) + Z - e) / 2 + e;%让所有Z都>e，Q=0对应的是nan
    Z_old = real(Z);
end
Z(~Q) = 0;

% % debug
s = sum(Z,'all');
if isnan(s) %|| ~isreal(s)
    error('Z有<=0的值')
end
end