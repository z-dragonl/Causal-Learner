
function [MB,test,time] = FBED_G2(Data,target,alpha,ns,p,maxK)
%
% FBED_G2 finds the Markov blanket of target node on discrete data
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


k_lim=1;    % Maximum Number of Runs

test=0;
S = [];

k_cur=0;

% Growing phase
while k_cur <= k_lim
    S_tmp = S;
    [S,test] = one_run(Data,target,S, alpha, ns, p, test);
    k_cur =k_cur+1;
    
    if length(S_tmp) == length(S)
        break;
    end
end

% Shrinking phase (same with IAMB)
tmp_MB = S;
for i=1:length(S)
    condset=mysetdiff(S,S(i));
    [pval]=my_g2_test(S(i),target,condset,Data,ns,alpha);
    test=test+1;
    if pval>alpha
        tmp_MB=mysetdiff(tmp_MB,S(i));
    end
end

S=tmp_MB;

MB=S;

time=toc(start);

end


function   [S, ntest]= one_run(Data,target,S, alpha, ns, p, ntest)

tmp_CanMB=[];
CanMB = mysetdiff([1:p],myunion(S,target));

while( ~isempty(CanMB) )
    
    tmp_dep = -10000;
    tmp_pval = 10000;
    
    for i=1:length(CanMB)      
        X = CanMB(i);
        
        ntest=ntest+1;
        [pval,dep]=my_g2_test(X,target,S,Data,ns,alpha);       
        if isnan(pval)
            tmp_CanMB = myunion(tmp_CanMB,X);
        else
            
            if pval>alpha              % early drop the nodes independent of T          
                tmp_CanMB = myunion(tmp_CanMB,X);  
            elseif (dep > tmp_dep)          
                
                tmp_dep = dep;
                tmp_pval = pval;
                Y = X;
                
            end
        end
    end
    
    if tmp_pval<=alpha
        S=[S Y];
        tmp_CanMB = myunion(tmp_CanMB,Y);
    else
        if tmp_pval>alpha 
            break;
        end
    end
    
    CanMB= mysetdiff(CanMB,tmp_CanMB);  
    
end

end

