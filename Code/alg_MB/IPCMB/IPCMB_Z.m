
function [MB,test,time] = IPCMB_Z(Data,target,alpha,samples,p,maxK)   
%
% IPCMB_Z finds the Markov blanket of target node on continuous data
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


test = 0;

% Recognize T' parents/children
[PC, ntest_PC, ~, Sepset] = PCsimple_Z(Data, target, alpha, samples, p, maxK);
test = test+ntest_PC;
MB = PC;

for i = 1:length(PC)
    % Recognize a true positive, and its 
    % parents/children as spouse candidates.
    
    X = PC(i);
    [CanSP, ntest_SP] = PCsimple_Z(Data, X, alpha, samples, p, maxK);
    test = test+ntest_SP;

    %---------------------------- Recognize a true positive
    if ~ismember(target,CanSP)
        MB = mysetdiff(MB,X);
        continue;
    end
    %-----------------------------end
    
    %---------------------------- Find spouse
    for j = 1:length(CanSP)
        Y = CanSP(j);  

        if ismember(Y,MB)||Y==target
            continue;
        end
        
        S =  myunion(Sepset{1,Y}, X);      % add X, then conditional independent => conditional dependent, so CI=0
        test = test +1;
        [CI]=my_fisherz_test(target,Y,S,Data,samples,alpha);       
        if isnan(CI)
            CI=0;
        end
        if(CI==0)    
            MB = myunion(MB,Y);  
        end     
    end
    %-----------------------------end
    
end

time=toc(start);


