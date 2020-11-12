

function [DAG,time,score] = MMHC_G2(Data,alpha,ns,p,maxK)
%
% MMHC_G2 learns a causal graph on discrete data
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
%       score is the score of the DAG found
%
%


if (nargin == 2)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end

start=tic;

score=0;

varNValues=ns;
sample=Data;



% create CIT-p-value estimator
CITPValueEstimator = gtestpvalueestimator(sample, varNValues);

% create CIT-RC applier
CITRCApplier = heuristicapplier(sample, varNValues);

% create d-separation determiner
CITDSepDeterminer = citdsepdeterminer(CITRCApplier, CITPValueEstimator);

% create and CA calculator
CITCACalculator = citcacalculator(CITRCApplier, CITPValueEstimator);



% create max-min local-to-global learner without symmetry correction
MMLGLearner = mmlglearner(...
    CITDSepDeterminer, ...
    CITCACalculator, ...
    'maxSepsetCard', 10, ...
    'symCorrEnabled', false);

% learn skeleton
skeleton = MMLGLearner.learnskeleton();


% create local scorer
LocalScorer = bdeulocalscorer(sample, varNValues);


% create candidate parent matrix
cpm = tril(skeleton + skeleton');


% create hill climber
HillClimber = hillclimber(LocalScorer, 'CandidateParentMatrix', cpm);


% Finally, we learn the structure of the network.

% learn structure
DAG = HillClimber.learnstructure();

time=toc(start);
