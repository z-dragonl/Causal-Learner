
function [DAG,time] = GES_G2(data,alpha,ns,p,maxK)
%
% GES_G2 learns a causal graph on discrete data
%
% INPUT :
%       Data is the data matrix
%       alpha is the significance level
%       ns is the size array of each node
%       p is the number of nodes
%       maxK is the maximum size of conditioning set
%
% OUTPUT:
%       DAG is the causal graph
%       time is the runtime of the algorithm
%
%


if (nargin == 2)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end


start=tic;

fprintf('IF the data is insurance_500.txt, \nplease remove column 16 of the data, which is constant.\n\n');
fprintf('For example, the input data needs to do this --> data(:,16)=[];\n');
fprintf('and the graph needs to do this --> graph(:,16)=[]; graph(16,:)=[];\n\n');

data=data-repmat(mean(data),size(data,1),1);
data=data*diag(1./std(data));

maxP = 5; % maximum number of parents when searching the graph
parameters.kfold = 10; % 10 fold cross validation
parameters.lambda = 0.01;
Record = GES_learn(data,1,maxP,parameters);
DAG=Record.G;

time=toc(start);

end