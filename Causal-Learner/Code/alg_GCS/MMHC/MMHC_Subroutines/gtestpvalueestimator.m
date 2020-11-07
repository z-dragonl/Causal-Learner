classdef(Sealed) gtestpvalueestimator < chi2citpvalueestimator
%GTESTPVALUEESTIMATOR G-test-p-value estimator.
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.CHI2CITPVALUEESTIMATOR is the
%   class of G-test-p-value estimators.

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

% References:
%   [1] http://en.wikipedia.org/wiki/G_test

methods

% constructor

function Obj = gtestpvalueestimator(varargin)
%GTESTPVALUEESTIMATOR Create G-test-p-value estimator.
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.GTESTPVALUEESTIMATOR/GTESTPVALUEESTIMATOR
%   creates a G-test-p-value estimator.
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.GTESTPVALUEESTIMATOR/GTESTPVALUEESTIMATOR
%   is called just like
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.CHI2CITPVALUEESTIMATOR/CHI2CITPVALUEESTIMATOR.
%
%   See also
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.CHI2CITPVALUEESTIMATOR/CHI2CITPVALUEESTIMATOR.

% call CHI2CITPVALUEESTIMATOR constructor
Obj = Obj@chi2citpvalueestimator(varargin{:});

end

end

methods(Access = protected)
    
% abstract method implementations

function stat = chi2stat(~, obs, exp)
    
% compute G statistic

terms = obs.*log(obs./exp);

% set NaN terms to zero
terms(isnan(terms)) = 0;
    
stat = 2*sum(terms);

end

end

end