classdef citpvalueestimator < handle
%CITPVALUEESTIMATOR Conditional-independence-test-p-value estimator.
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CITPVALUEESTIMATOR is the class of
%   conditional-independence-test-p-value estimators. A
%   conditional-independence-test-p-value estimator estimates p-values of
%   testing conditional-independence relationships among a set of
%   variables.

% Copyright 2010-2012 Mens X Machina
% 
% This file is part of Mens X Machina Common Toolbox.
% 
% Mens X Machina Common Toolbox is free software: you can redistribute it
% and/or modify it under the terms of the GNU General Public License
% alished by the Free Software Foundation, either version 3 of the License,
% or (at your option) any later version.
% 
% Mens X Machina Common Toolbox is distributed in the hope that it will be
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along
% with Mens X Machina Common Toolbox. If not, see
% <http://www.gnu.org/licenses/>.

properties(Abstract, SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer
    
end

methods(Abstract)
    
%CITPVALUE Estimate p-value of hypothesis test of conditional independence.
%   P = CITPVALUE(CITPESTOBJ, I, J, K), where CITPESTOBJ is a
%   conditional-independence-test-p-value estimator, I and J is the linear
%   index of CITPESTOBJ variable X and Y, respectivelly, and K are the
%   linear indices of set Z of CITPESTOBJ variables, estimates the p-value
%   of testing the hypothesis that X and Y are independent given Z. I and J
%   are numeric integers in range [1, M], where M is the number of
%   CITPESTOBJ variables. K is a numeric row vector of integers in range
%   [1, M]. P is a numeric real scalar in [0, 1] or NaN if, for some
%   reason, the p-value cannot be estimated.
%
%   [P T] = CITPVALUE(CITPESTOBJ, I, J, K) also returns the statistic T of
%   the test. T is a numeric real scalar or NaN if, for some reason, the
%   p-value cannot be estimated.
%
%   Example:
%
%       % pValueEstimator is a conditional-independence-test-p-value
%       % estimator
%
%       % estimate the p-value of testing Ind(1, 2| [3 4])
%       p = citpvalue(pValueEstimator, 1, 2, [3 4]);
[p t] = citpvalue(Obj, i, j, k);

end

methods(Access = protected)
    
% input parsers    

function parsecitpvalueinput(Obj, i, j, k)
%PARSECITPVALUEINPUT Parse ORG.MENSXMACHINA.STATS.TESTS.CI.CITPVALUEESTIMATOR/CITPVALUE input.
%   PARSECITPVALUEINPUT(CITPESTOBJ, ...), when CITPESTOBJ is a
%   conditional-independence-test-p-value estimator, throws an error if its
%   input is not valid input for
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CITPVALUEESTIMATOR/CITPVALUE.
%
%   See also ORG.MENSXMACHINA.STATS.TESTS.CI.CITPVALUEESTIMATOR/CITPVALUE.
    
    validateattributes(i, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    validateattributes(j, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    validateattributes(k, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    
end

end

end