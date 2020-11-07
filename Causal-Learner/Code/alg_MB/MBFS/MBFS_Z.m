
function [MB,test,time] = MBFS_Z(Data,target,alpha,samples,p,maxK)
%
% MBFS_Z finds the Markov blanket of target node on continuous data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       samples is the number of data samples
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
   [samples,p]=size(Data);
   maxK=3;
end

start=tic;


test=0;

[pc,ntest1]=FCBF_PC_Z(Data,target,alpha, samples, p, maxK);
test=test+ntest1;

MB=pc;

for i=1:length(pc)
    [pc_tmp,ntest2]=FCBF_PC_Z(Data,pc(i),alpha, samples, p, maxK);
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
                CI=my_fisherz_test(pc_tmp(j),target,Z,Data,samples,alpha);
                
                if CI
                    
                    if ~ismember(pc(i),Z)
                        test=test+1;
                        CI=my_fisherz_test(pc_tmp(j),target,myunion(Z,pc(i)),Data,samples,alpha);
                        
                        if isnan(CI)||CI==0
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





