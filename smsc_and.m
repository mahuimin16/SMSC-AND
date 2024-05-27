function [obj,resiter,res] = smsc_and(X,Y,lambda,cri)
% X      : n*di
maxIter = 20 ; % the number of iterations
numclass = length(unique(Y)) ;
numview = length(X);
numsample = size(Y,1);
D = cell(numview,1);  
Z =  cell(numview,1); 
Q = Z;  % n * n
Beta = zeros(numview,1);        % v * 1
Qc =  ones(numsample,numsample);
C =  zeros(numsample,numsample);
miu =  cell(numview,1); 
niu =  cell(numview,1); 
miuC =  []; 
niuC =  []; 
%% initialize
for i = 1:numview
    di = size(X{i},2);
    X{i} = mapminmax(X{i}',0,1); % turn into d*n
    Beta(i) = 1/numview;
    D{i} = pdist(X{i}').^2;
    D{i} = squareform(D{i});
    Q{i} = ones(numsample,numsample);
    Q{i} = ones(numsample,numsample)-diag(diag(Q{i}));
    Z{i} = 1/(numsample-1).*ones(numsample,numsample);
    Z{i} = Z{i} - diag(diag(Z{i}));
end
C(:) = 1/(numsample-1);
C = C - diag(diag(C));
Qc = Qc - diag(diag(Qc));


%% Iteriation
flag = 1;
iter = 0;
obj = [];
while flag
    iter = iter +1;
    term1 = 0;
    term2 = 0;
    term3 = 0;
    term4 = 0;
    for i =1: numview
        term1 = term1 + sum(D{i} .* Z{i},'all' ) ;
        Z_temp = Z{i};
        Z_temp (Z{i}<=1/numsample/100) = 1;
        lnZ = Z_temp.* log(Z_temp) ;
        term2 =term2 + lambda * sum(lnZ,'all') ;
        term3 =  term3 + Beta(i)^2 .* sum((C - Z{i}).^2,'all') ;
    end
    C_temp = C;
    C_temp (C<=1/numsample/100) = 1;
    lnC = C_temp.* log(C_temp) ;
    term4 = term4 + lambda * sum(lnC,'all');
    term = term1 + term2 + term3 + term4;
obj(iter,:) = [term term1 term2 term3 term4];
    if (iter>5)&&(abs(obj(iter-1)-obj(iter))/abs(obj(3)-obj(4))<1e-2 || iter>maxIter || abs(obj(iter)) < 1e-10)
        flag = 0;
    end
    
%% Update Z^{(p)}

    for i=1:numview
        [Z{i},miu{i},niu{i}] = updateZ(D{i}-2*Beta(i)^2.*C,lambda,Beta(i)^2,Q{i},Z{i},miu{i},niu{i});
    end


%% Update Q^{(p)}
if iter == 1
    for i=1:numview
        [~,~,~,Q{i}] = criterion(Z{i},cri) ;
        Q{i} = Q{i} + Q{i}' > 0;
    end
end

%% Update C
conZ = zeros(numsample,numsample);
    for iv = 1:numview
        conZ = conZ + Beta(iv)^2.*Z{iv};
    end
[C,miuC,niuC] = updateZ(-2*conZ,lambda,Beta'*Beta,Qc,C,miuC,niuC);
%% Update Qc
if iter == 1
        [~,~,~,Qc] = criterion(C,cri) ;
        Qc = Qc + Qc' > 0;
end
%% Update Beta
M = zeros(numview,1);
    for iv = 1:numview
        M(iv) = sum((C-Z{iv}).^2,'all'); %+ gamma*sum( lnZ{iv}, 'all') + sum(min(Cri{iv}))
    end
    M_inv = 1 ./ M;
    Beta = M_inv / (sum(M_inv));

%%
    [label,res] = SpectralClustering(C.*Qc,numclass,Y);
    resiter(iter,:) = res;
end
    

 