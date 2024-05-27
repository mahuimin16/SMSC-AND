%--------------------------------------------------------------------------
% This function takes an adjacency matrix of a graph and computes the 
% clustering of the nodes using the spectral clustering algorithm of 
% Ng, Jordan and Weiss.
% CMat: NxN adjacency matrix
% n: number of groups for clustering
% groups: N-dimensional vector containing the memberships of the N points 
% to the n groups obtained by spectral clustering
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function [label,maxres] = SpectralClustering(CKSym,n,gt)

warning off;
N = size(CKSym,1);
 MAXiter = 100; % Maximum number of iterations for KMeans 
 REPlic = 20; % Number of replications for KMeans
maxiter = 50;

% Normalized spectral clustering according to Ng & Jordan & Weiss
% using Normalized Symmetric Laplacian L = I - D^{-1/2} W D^{-1/2}
DN = diag( 1./sqrt(sum(CKSym)) );
LapN = speye(N) - DN * CKSym * DN;
[uN,sN,vN] = svd(LapN);
kerN = vN(:,N-n+1:N);
for i = 1:N
    kerNS(i,:) = kerN(i,:) ./ norm(kerN(i,:)+eps);
end
% for i = 1:maxiter
% l = kmeans(kerNS,n,'MaxIter',100, 'Replicates',1);
%  [ress(i,:)] = Clustering8Measure(gt,l);
% end
% result = max(ress,[],1);
for rep= 1: maxiter
% pY = kmeans(kerNS,n,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
  pY = kmeans(kerNS,n);
 res(rep, : ) = Clustering8Measure(gt, pY);
 plabel(: , rep) = pY ;  
end
[maxres,I] = max(res(1:maxiter, : ));
label = plabel(:,I);
