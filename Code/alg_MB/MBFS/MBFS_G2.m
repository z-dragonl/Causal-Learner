
function [MB,test,time] = MBFS_G2(Data,target,alpha,ns,p,maxK)
%
% MBFS_G2 finds the Markov blanket of target node on discrete data
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

[pc,ntest1]=FCBF_PC_G2(Data,target,alpha, ns, p, maxK);
test=test+ntest1;

MB=pc;

for i=1:length(pc)
    [pc_tmp,ntest2]=FCBF_PC_G2(Data,pc(i),alpha, ns, p, maxK);
    test=test+ntest2;

    for j=1:length(pc_tmp)
        
        if ismember(pc_tmp(j),pc)||pc_tmp(j)==target
            continue;
        end
                
        CanPC=pc;
        cutSetSize = 0;
        
        break_flag=0;

        while length(CanPC) >= cutSetSize &&cutSetSize<=maxK
            
            SS = subsets1(CanPC, cutSetSize); 
            for si=1:length(SS)
                Z = SS{si};

                test=test+1;
                pval=my_g2_test(pc_tmp(j),target,Z,Data,ns,alpha);
                
                if pval>alpha
                    
                    if ~ismember(pc(i),Z)
                        test=test+1;
                        pval=my_g2_test(pc_tmp(j),target,myunion(Z,pc(i)),Data,ns,alpha);
                        
                        if isnan(pval)||pval<=alpha
                            MB=myunion(MB,pc_tmp(j));
                        end
                    end
                    
                    break_flag=1;
                    break;
                end
            end
            if break_flag
                break;
            end
            cutSetSize = cutSetSize + 1;
        end
        
    end
end

time=toc(start);





