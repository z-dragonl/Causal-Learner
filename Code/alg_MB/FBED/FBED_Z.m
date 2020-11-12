
function [MB,test,time] = FBED_Z(Data,target,alpha,samples,p,maxK)
%
% FBED_Z finds the Markov blanket of target node on continuous data
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


k_lim=1;    % Maximum Number of Runs

test=0;
S = [];

k_cur=0;

% Growing phase
while k_cur <= k_lim
    S_tmp = S;
    [S,test] = one_run(Data,target,S, alpha, samples, p, test);
    k_cur =k_cur+1;
    
    if length(S_tmp) == length(S)
        break;
    end
end

% Shrinking phase (same with IAMB)
tmp_MB = S;
for i=1:length(S)
    condset=mysetdiff(S,S(i));
    [CI]=my_fisherz_test(S(i),target,condset,Data,samples,alpha);
    test=test+1;
    if CI
        tmp_MB=mysetdiff(tmp_MB,S(i));
    end
end

S=tmp_MB;

MB=S;

time=toc(start);

end


function   [S, ntest]= one_run(Data,target,S, alpha, samples, p, ntest)

tmp_CanMB=[];
CanMB = mysetdiff([1:p],myunion(S,target));

while( ~isempty(CanMB) )
    
    tmp_dep = -10000;
    tmp_CI = 10000;
    
    for i=1:length(CanMB)      
        X = CanMB(i);
        
        ntest=ntest+1;
        [CI,dep]=my_fisherz_test(X,target,S,Data,samples,alpha);        
        if isnan(CI)
            tmp_CanMB = myunion(tmp_CanMB,X);
        else
            
            if CI              % early drop the nodes independent of T          
                tmp_CanMB = myunion(tmp_CanMB,X);  
            elseif (dep > tmp_dep)         
                
                tmp_dep = dep;
                tmp_CI = CI;
                Y = X;
                
            end
        end
    end
    
    if tmp_CI==0
        S=[S Y];
        tmp_CanMB = myunion(tmp_CanMB,Y);
    else
        if tmp_CI
            break;
        end
    end
    
    CanMB= mysetdiff(CanMB,tmp_CanMB);  
    
end

end

