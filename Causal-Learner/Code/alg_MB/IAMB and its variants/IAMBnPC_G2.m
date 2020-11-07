
function [MB,test,time] = IAMBnPC_G2(Data,target,alpha,ns,p,maxK)
%
% IAMBnPC_G2 finds the Markov blanket of target node on discrete data
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


MB=[];

test=0;

tmp_CanMB=[];
CanMB = mysetdiff([1:p],target);

%---------------------------------------------------------
%  add true positives to MB

while( ~isempty(CanMB) )
    
    tmp_dep = -10000;
    tmp_pval = 10000;
    
    for i=1:length(CanMB)      
        X = CanMB(i);
        
        test=test+1;
        [pval,dep]=my_g2_test(X,target,MB,Data,ns,alpha);      
        
        if isnan(pval)
            tmp_CanMB = [tmp_CanMB,X];
        else
            if (dep > tmp_dep)          
                
                tmp_dep = dep;
                tmp_pval = pval;
                Y = X;
                
            end
        end
    end
    
    if tmp_pval<=alpha
        MB=[MB Y];
        tmp_CanMB=[tmp_CanMB,Y];
        CanMB= mysetdiff(CanMB,tmp_CanMB);
    else
        
        if tmp_pval>alpha 
            break;
        end
    end
end


% -----------------------------------------------------------
% remove false positives from MB

tmp_MB = MB;

for i=length(tmp_MB):-1:1
    Y=tmp_MB(i);
    SZ = mysetdiff(MB,Y);
    cutSetSize = 0;
    
    break_flag=0;
    while length(SZ) >= cutSetSize&&cutSetSize<=maxK
        SS = subsets1(SZ, cutSetSize);   
        for si=1:length(SS)
            S = SS{si};
            test=test+1;
            [pval]=my_g2_test(Y,target,S,Data,ns,alpha);       
            if isnan(pval)        
                CI=0;
            else
                if pval<=alpha
                    CI=0;
                else
                    CI=1;
                end
            end
            if(CI==1)             
                MB=SZ;
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

MB = sort(MB);

time=toc(start);





