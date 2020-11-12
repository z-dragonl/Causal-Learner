classdef(Sealed) bdeulocalscorer < localscorer
%BDEULOCALSCORER BDeu local scorer.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.BDEU.BDEULOCALSCORER
%   is the class of BDeu local scorers.

% Copyright 2010-2012 Mens X Machina
% 
% This file is part of Mens X Machina Probabilistic Graphical Model
% Toolbox.
% 
% Mens X Machina Probabilistic Graphical Model Toolbox is free software:
% you can redistribute it and/or modify it under the terms of the GNU
% General Public License alished by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version.
% 
% Mens X Machina Probabilistic Graphical Model Toolbox is distributed in
% the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
% the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE. See the GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with Mens X Machina Probabilistic Graphical Model Toolbox. If not, see
% <http://www.gnu.org/licenses/>.

% References:
% [1] Heckerman, D. E., Geiger, D., & Chickering, D. M. (1995). Learning
%     Bayesian networks: The combination of knowledge and statistical 
%     sample. Machine Learning, 20, 197-243.

properties(SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer

end

properties(GetAccess = private, SetAccess = immutable)
    
sample % sample -- an M-by-N numeric matrix of positive integers

varNValues % variable numbers of values -- an 1-by-N numeric vector of positive integers

equivalentSampleSize % equivalent sample size -- a numeric nonnegative integer
    
end

methods

% constructor

function Obj = bdeulocalscorer(sample, varNValues, equivalentSampleSize)
%BDEULOCALSCORER Create BDeu local scorer.
%   OBJ =
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.BDEU.BDEULOCALSCORER(D,
%   VARNVALUES) creates a BDeu local scorer with sample D and variable
%   numbers of values VARNVALUES. D is an M-by-N numeric matrix of positive
%   integers, where M is the number of observations and N is the number of
%   variables in the sample. Each element of D is the linear index of the
%   value of the corresponding variable for the corresponding observation.
%   VARNVALUES is an 1-by-N numeric vector of positive integers. Each
%   element of VARNVALUES is the number of values of the corresponding
%   variable.
%
%   OBJ =
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.BDEU.BDEULOCALSCORER(D,
%   VARNVALUES, ESS) uses equivalent sample size ESS. ESS is a is a numeric
%   nonnegative integer.

    % parse input
    validateattributes(sample, {'numeric'}, {'2d', 'positive', 'integer'});
    validateattributes(varNValues, {'numeric'}, {'size', [1 size(sample, 2)], 'positive', 'integer'});    
    
    if nargin < 3
        equivalentSampleSize = 10;
    else
        validateattributes(equivalentSampleSize, {'numeric'}, {'scalar', 'nonnegative', 'integer'});
    end
    
    % set properties
    Obj.nVars = size(sample, 2);
    Obj.sample = sample;
    Obj.varNValues = varNValues;
    Obj.equivalentSampleSize = equivalentSampleSize;

end

end

methods

% abstract method implementations

function score = loglocalscore(Obj, i, j)

    
    % select PA(X) number of values
    xPANValues = Obj.varNValues(j);
    
    % select X number of values
    xNValues = Obj.varNValues(i);
    
    % select FA(X) sample
    xFASample = Obj.sample(:, [j i]);
    
    sampleSize = size(xFASample, 1);
    
    % compute counts
    N = accumarray(xFASample, ones(1, sampleSize), makesize([xPANValues xNValues]));

    nXPAValues = prod(xPANValues);

    N_ijk = reshape(N, nXPAValues, xNValues);
    N_ij = sum(N_ijk, 2);
    N_prime_ijk = Obj.equivalentSampleSize*ones(nXPAValues, xNValues)./(nXPAValues*xNValues);
    N_prime_ij = sum(N_prime_ijk, 2);

    % compute log prior [1] / # nodes
    delta = length([j i]) - 1; % because the prior network is the empty network
    kappa = 1/(Obj.equivalentSampleSize + 1);
    prior = delta*log(kappa);

    % compute log likelihood
    likelihood = sum(gammaln(N_prime_ij) - gammaln(N_prime_ij + N_ij)) + ...
                 sum(sum(gammaln(N_prime_ijk + N_ijk) - gammaln(N_prime_ijk)));
             
    % add them         
    score = prior + likelihood;

end

end

end