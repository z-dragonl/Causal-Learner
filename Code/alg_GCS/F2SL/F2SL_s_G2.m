

function [DAG,time] = F2SL_s_G2(Data,alpha,ns,p,maxK)
%
% F2SL_s_G2 learns a causal graph on discrete data
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

varNValues=ns;
sample=Data;

skeleton=zeros(p,p);
for i=1:p
    PC=FCBF_PC_G2(sample,i,alpha,ns,p,maxK);
    skeleton(i,PC)=1;
end

skeleton=sparse(skeleton);

% create local scorer
LocalScorer = bdeulocalscorer(sample, varNValues);


% create candidate parent matrix
cpm = tril(sign(skeleton + skeleton'));


% create hill climber
HillClimber = hillclimber(LocalScorer, 'CandidateParentMatrix', cpm);


% Finally, we learn the structure of the network.

% learn structure
DAG = HillClimber.learnstructure();

% draw_graph(DAG);

time=toc(start);
