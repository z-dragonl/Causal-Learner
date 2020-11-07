
function [DAG,time] = GES_Z(data,alpha,samples,p,maxK)
%
% GES_Z learns a causal graph on continuous data
%
% INPUT :
%       Data is the data matrix
%       alpha is the significance level
%       samples is the number of data samples
%       p is the number of nodes
%       maxK is the maximum size of conditioning set
%
% OUTPUT:
%       DAG is the causal graph
%       time is the runtime of the algorithm
%
%


if (nargin == 2)
   [samples,p]=size(Data);
   maxK=3;
end


start=tic;

data=data-repmat(mean(data),size(data,1),1);
data=data*diag(1./std(data));

maxP = 5; % maximum number of parents when searching the graph
parameters.kfold = 10; % 10 fold cross validation
parameters.lambda = 0.01;
Record = GES_learn(data,1,maxP,parameters);
DAG=Record.G;

time=toc(start);

end