
function [MB,test,time] = GS_G2(Data,target,alpha,ns,p,maxK)
%
% GS_G2 finds the Markov blanket of target node on discrete data
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

CanMB = mysetdiff([1:p],target);

%---------------------------------------------------------
%  add true positives to MB


while( ~isempty(CanMB) )
    
    break_flag=1;
    
    for i=1:length(CanMB)      
        X = CanMB(i);
        test=test+1;
        [pval]=my_g2_test(X,target,MB,Data,ns,alpha);       

        if isnan(pval)
            CI=1;
        else
            if pval<=alpha
                CI=0;
            else
                CI=1;
            end
        end
        
        if CI==0
            break_flag=0;
            MB = [MB,X];
        end
    end
    
    if break_flag   % all nodes in CanMB are independent of T
        break;
    end
    
    CanMB = mysetdiff(CanMB,MB);

end

% -----------------------------------------------------------
% remove false positives from MB

tmp_MB = MB;
for i=1:length(MB)
    condset=mysetdiff(MB,MB(i));
    [pval]=my_g2_test(MB(i),target,condset,Data,ns,alpha);
    test=test+1;
    if pval>alpha
        tmp_MB=mysetdiff(tmp_MB,MB(i));
    end
end

MB = tmp_MB;

time=toc(start);





