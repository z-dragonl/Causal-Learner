
function [MB,test,time] = interIAMBnPC_Z(Data,target,alpha,samples,p,maxK)
%
% interIAMBnPC_Z finds the Markov blanket of target node on continuous data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       samples is the size array of each node
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



while( ~isempty(CanMB) )
    
    %---------------------------------------------------------
    %  add true positives to MB
    
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
    
    % -----------------------------------------------------------
    % remove false positives from MB
    
    
    tmp_MB = MB;
    
    last_break_flag=0;
    
    for i=length(tmp_MB):-1:1

        SZ = mysetdiff(MB,tmp_MB(i));
        cutSetSize = 0;
        
        break_flag=0;
        while length(SZ) >= cutSetSize&&cutSetSize<=maxK
            SS = subsets1(SZ, cutSetSize);   
            for si=1:length(SS)
                S = SS{si};
                
                if i~=length(tmp_MB)&&~ismember(tmp_MB(i),S)
                    continue;
                end
                
                test=test+1;
                [CI]=my_fisherz_test(tmp_MB(i),target,S,Data,samples,alpha);      
                if isnan(CI)       
                    CI=0;
                end
                if(CI==1)         
                    MB=SZ;
                    break_flag=1;
                    
                    if i==length(tmp_MB)
                        last_break_flag=1;
                    end
                    
                    break;
                end
            end
            if break_flag||last_break_flag
                break;
            end
            cutSetSize = cutSetSize + 1;
        end
        
        if last_break_flag
            break;
        end
    end
    
end


MB = sort(MB);

time=toc(start);





