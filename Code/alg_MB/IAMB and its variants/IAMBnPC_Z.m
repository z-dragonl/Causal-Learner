
function [MB,test,time] = IAMBnPC_Z(Data,target,alpha,samples,p,maxK)
%
% IAMBnPC_Z finds the Markov blanket of target node on continuous data
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


MB=[];

test=0;

tmp_CanMB=[];
CanMB = mysetdiff([1:p],target);

%---------------------------------------------------------
%  add true positives to MB

while( ~isempty(CanMB) )
    
    tmp_dep = -10000;
    tmp_CI = 10000;
    
    for i=1:length(CanMB)      
        X = CanMB(i);
        
        test=test+1;
        [CI,dep]=my_fisherz_test(X,target,MB,Data,samples,alpha);      
        
        if isnan(CI)
            tmp_CanMB = [tmp_CanMB,X];
        else
            if (dep > tmp_dep)         
                
                tmp_dep = dep;
                tmp_CI = CI;
                Y = X;
                
            end
        end
    end
    
    if tmp_CI==0
        MB=[MB Y];
        tmp_CanMB=[tmp_CanMB,Y];
        CanMB= mysetdiff(CanMB,tmp_CanMB);
    else
        
        if tmp_CI
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
            [CI]=my_fisherz_test(Y,target,S,Data,samples,alpha);       
            if isnan(CI)       
                CI=0;
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





