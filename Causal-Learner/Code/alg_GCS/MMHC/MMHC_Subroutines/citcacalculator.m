classdef(Sealed) citcacalculator < cacalculator
%CITCACALCULATOR Conditional-independence-test-based conditional-association calculator.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CIT.CITCACALCULATOR is the
%   class of conditional-independence-test-based conditional-association
%   calculators. A conditional-independence-test-based
%   conditional-association calculator performs a hypothesis test of
%   conditional independence for each conditional association to be
%   calculated. If the test is not completed, both associations are set to
%   Inf. If the test is completed and the null hypothesis is rejected, the
%   primary association is set to 1/P, where P is the p-value of the test,
%   and the secondary association is set to ABS(T), where T is the
%   statistic of the test. If the test is completed and the null hypothesis
%   cannot be rejected, both associations are set to 0.
%
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CIT.CITCACALCULATOR/CA
%   triggers a citPerformed event with
%   conditional-independence-test-performed data when a test is performed.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CACALCULATOR,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.CIT.CITPERFORMEDDATA.

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
% [1] I. Tsamardinos, L.E. Brown, C.F. Aliferis. "The Max-Min Hill-Climbing
%     Bayesian Network structure Learning Algorithm", Machine Learning
%     Journal; 65: 31-78

properties(SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer
    
CITPValueEstimator % conditional-independence-test-p-value estimator
CITRCApplier % conditional-independence-test-reliability-criterion applier

type1ErrorRate % type I error rate -- a numeric real scalar in range [0, 1]
    
end

events

citPerformed % conditional-independence test performed

end

methods

% constructor

function Obj = citcacalculator(CITRCApplier, CITPValueEstimator, alpha)
%CITCACALCULATOR Create conditional-independence-test-based conditional-association calculator.
%   CITPERFOBJ =
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CIT.CITCACALCULATOR(CITRCAPPOBJ,
%   CITPESTOBJ) creates a
%	conditional-independence-test-based conditional-association calculator
%	with reliability-criterion applier CITRCAPPOBJ, p-value estimator
%	CITPESTOBJ, and type I error rate 0.05. CITPESTOBJ and CITRCAPPOBJ have
%	the same number of variables.
%
%   CITPERFOBJ =
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CIT.CITCACALCULATOR(CITRCAPPOBJ,
%   CITPESTOBJ, ALPHA) creates a calculator with type I error rate ALPHA.
%   ALPHA is a numeric real scalar in range [0, 1].

    % parse input
    assert(isa(CITRCApplier, 'citrcapplier') && isscalar(CITRCApplier));
    assert(isa(CITPValueEstimator, 'citpvalueestimator') && isscalar(CITPValueEstimator) && CITPValueEstimator.nVars == CITRCApplier.nVars);
    
    if nargin < 3
        alpha = 0.05;
    else
        validateattributes(alpha, {'numeric'}, {'real', 'scalar', '>=', 0, '<=', 1});
    end
    
    % set properties
    
    Obj.nVars = CITPValueEstimator.nVars;
    
    Obj.CITRCApplier = CITRCApplier;
    Obj.CITPValueEstimator = CITPValueEstimator;
    
    Obj.type1ErrorRate = alpha;

end

% abstract method implementations

function [ca ca2] = ca(Obj, i, j, k)
    
    import citperformeddata;

    if isreliablecit(Obj.CITRCApplier, i, j, k) % test is reliable, attempt it

        % calculate test p-value
        [p stat] = citpvalue(Obj.CITPValueEstimator, i, j, k);

        % compare p-value with alpha

        if p > Obj.type1ErrorRate % independent

            ca = 0;
            ca2 = 0;

        else % dependent

            ca = 1/p;
            ca2 = abs(stat);

        end

    else
        
        p = NaN;
        stat = NaN;

        if isempty(k)
            ca = 0;
            ca2 = 0;
        else
            ca = Inf;
            ca2 = Inf;
        end

    end
            
    notify(Obj, 'citPerformed', citperformeddata(i, j, k, p, stat));

end

end

end