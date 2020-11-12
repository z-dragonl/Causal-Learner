
function [MB,test,time] = HITONMB_G2(Data,target,alpha,ns,p,maxK)
%
% HITONMB_G2 finds the Markov blanket of target node on discrete data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       ns is the size array of each node
%       p is the number of nodes
%       maxK is the maximum size of conditioning set
%
% OUTPUT:
%       MB is the Markov blanket of the target
%       test is the number of conditional independence tests
%       time is the runtime of the algorithm
%
%


if (nargin == 3)
    ns=max(Data);
    [~,p]=size(Data);
    maxK=3;
end

start=tic;

test=0;

sp=[];
MB=[];


[pc,ntest1,~,sepset]=HITONPC_G2(Data,target,alpha, ns, p, maxK);
test=test+ntest1;
MB=[MB pc];


for i=1:length(pc)
    
    [pc_tmp,ntest2]=HITONPC_G2(Data,pc(i),alpha, ns, p, maxK);
    test=test+ntest2;
    for j=1:length(pc_tmp)
        
        if isempty(find(pc==pc_tmp(j), 1))&& pc_tmp(j)~=target && isempty(find(sepset{pc_tmp(j)}==pc(i), 1))
            
            [pval]=my_g2_test(pc_tmp(j),target,[sepset{pc_tmp(j)},pc(i)],Data,ns,alpha);
            
            if isnan(pval)
                CI=0;
            else
                if pval<=alpha
                    CI=0;
                else
                    CI=1;
                end
            end
            
            test=test+1;
            if CI==0
                sp=myunion(sp,pc_tmp(j));
            end
        end
    end
end

MB=myunion(MB,sp);
time=toc(start);