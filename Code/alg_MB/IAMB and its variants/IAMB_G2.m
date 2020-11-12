
function [MB,test,time] = IAMB_G2(Data,target,alpha,ns,p,maxK)
%
% IAMB_G2 finds the Markov blanket of target node on discrete data
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
    else
        if tmp_pval>alpha 
            break;
        end
    end
    
    CanMB= mysetdiff(CanMB,tmp_CanMB);
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

MB = sort(tmp_MB);

time=toc(start);





