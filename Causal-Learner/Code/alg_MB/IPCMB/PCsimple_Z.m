
function [PC,test,time,sepset] = PCsimple_Z(Data,target,alpha,samples,p,maxK)      
%
% PCsimple_Z finds the parents and children of target node on continuous data
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
%       PC is the parents and children of the target
%       test is the number of conditional independence tests
%       time is the runtime of the algorithm
%       sepset: the condition sets for each node to make it independt to the target
%
%


if (nargin == 3)
   [samples,p]=size(Data);
   maxK=3;
end

start=tic;

test=0;


sepset = cell(1,p);

ADJT = mysetdiff(1:p,target);
ADJT_length = length(ADJT);
cutSetSize = 0;
tmp_ADJT = ADJT;

while ADJT_length > cutSetSize&&cutSetSize<=maxK
    for i=1:length(ADJT)
        X = ADJT(i);
        nbrs = mysetdiff(tmp_ADJT, X);
        
        SS = subsets1(nbrs, cutSetSize);   
        for si=1:length(SS)
            S = SS{si};
            test=test+1;
            [CI]=my_fisherz_test(X,target,S,Data,samples,alpha);        
            if isnan(CI)
                CI=0;
            end            
            
            if(CI==1)           
                tmp_ADJT = mysetdiff(tmp_ADJT,X);
                ADJT_length = ADJT_length-1;
                sepset{1,X} = S;
                break;          
            end
        end
        
    end

    ADJT=tmp_ADJT;
    
    cutSetSize = cutSetSize + 1;
end

PC=ADJT;

time=toc(start);


