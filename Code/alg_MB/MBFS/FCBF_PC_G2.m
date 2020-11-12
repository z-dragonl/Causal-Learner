
function [PC,test,time] = FCBF_PC_G2(Data,target,alpha,ns,p,maxK)    
%
% FCBF_PC_G2 finds the parents and children of target node on discrete data
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
%       PC is the parents and children of the target
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

threshold = 0.05; % threshod of the feature selection method FCBF

test=0;

train_data=Data(:,mysetdiff(1:p,target));
train_label=Data(:,target);


% Map=unique([1:target-1,target+1:p]);


if target~=p 
    Map=unique([1:target-1,target+1:p]);
else
    Map=1:target;
end

[FCBF_PC,ntest_PC]=FCBF(train_data,train_label,threshold);

test=test+ntest_PC;

PC=sort(Map(FCBF_PC));

time=toc(start);





