clear;
% clc;
warning off;
addpath(genpath('./'));

%% dataset

ds = {'MSRCV1'};
% dsPath = '.\datasets\';
dsPath = 'D:\Master\Research\Code\MVC-E\datasets\';
resultdir = '.\res\';
metric = {'ACC','nmi','Purity','Fscore','Precision','Recall','AR','Entropy'};
lambda = [8];
cri = ["in_e"];

% cri = ["in_e","diff"];


for dsi =1:length(ds)
    dataName = ds{dsi}; disp(dataName);
    load(strcat(dsPath,dataName));
%     X= data';
%     Y= truth ;
    k = length( unique(Y));
    n = length(Y);
    %%
    for dcri = 1:length(cri)
    for id = 1:length(lambda)      
        tic;
        [obj,resiter,res] = smsc_and(X,Y,lambda(id),cri(dcri)); 
        timer(id,dcri)  = toc;
        resall{id,dcri} = res;
        resitall{id,dcri} = resiter;
        objall{id,dcri} = obj;
    end
    end
    save([resultdir, char(dataName),'_result.mat'], 'resall', 'resitall','objall','timer');
end


