 function [Cri,Z_bestloca,Z_best,Q] = criterion(Z,cri)
 if nargin < 2
     cri = ["in_e"];
end
[S,I]=sort(Z,'descend');
n = size(S,2);
Q = zeros(n,n);
Cri = zeros(n,n);
if cri == "in_e"
    for i=1:n
%         kn = sum(S(:,i)>0);
            kn = n-1;
            for k=1:kn
            Cri(k,i)=-k*(kn-k)/kn/kn*(1/k*sum(S(1:k,i))-1/(kn-k)*sum(S((k+1):kn,i))).^2;
            end
        [~,Z_bestloca(i)]=min(Cri(:,i));
        Z_best(i)=S(Z_bestloca(i),i);
        Q(1:n,i) =(Z(1:n,i)>=Z_best(i));
    end
end
if cri == "diff"
    for i=1:n
        for k=1:n-2
            Cri(k,i)=S(k+1)-S(k,i);
        end
        [~,Z_bestloca(i)]=min(Cri(:,i));
        Z_best(i)=S(Z_bestloca(i),i);
        Q(1:n,i) =(Z(1:n,i)>=Z_best(i));
    end
end