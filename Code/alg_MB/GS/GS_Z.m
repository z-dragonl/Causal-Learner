
function [MB,test,time] = GS_Z(Data,target,alpha,samples,p,maxK)
%
% GS_Z finds the Markov blanket of target node on continuous data
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

CanMB = mysetdiff([1:p],target);

%---------------------------------------------------------
%  add true positives to MB


while( ~isempty(CanMB) )
    
    break_flag=1;
    
    for i=1:length(CanMB)      
        X = CanMB(i);
        test=test+1;
        [CI]=my_fisherz_test(X,target,MB,Data,samples,alpha);      

        if isnan(CI)
            CI=1;
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
    [CI]=my_fisherz_test(MB(i),target,condset,Data,samples,alpha);
    test=test+1;
    if CI
        tmp_MB=mysetdiff(tmp_MB,MB(i));
    end
end

MB = tmp_MB;

time=toc(start);





