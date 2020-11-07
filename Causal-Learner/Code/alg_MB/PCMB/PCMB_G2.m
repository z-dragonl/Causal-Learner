
function [MB,test,time] = PCMB_G2(Data,target,alpha,ns,p,maxK)
%
% PCMB_G2 finds the Markov blanket of target node on discrete data
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

already_calculated_PCD=ones(1,p);
all_PCD=cell(1,p);


[pc,ntest1,~,sepset]=GetPC_G2(Data,target,alpha, already_calculated_PCD, all_PCD, ns, p, maxK);
% [pc,ntest1,~,sepset]=GetPCD_G2(Data,target,alpha, ns, p, maxK);
 test=test+ntest1;
 MB=[MB pc];

 
for i=1:length(pc)

    [pc_tmp,ntest2]=GetPC_G2(Data,pc(i),alpha, already_calculated_PCD, all_PCD, ns, p, maxK);
%     [pc_tmp,ntest2]=GetPCD_G2(Data,pc(i),alpha, ns, p, maxK);
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