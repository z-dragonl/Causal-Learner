classdef(Sealed) citdsepdeterminer < dsepdeterminer
%CITDSEPDETERMINER Conditional-independence-test-based d-separation determiner.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.CIT.CITDSEPDETERMINER is the class
%   of conditional-independence-test-based d-separation determiners. A
%   conditional-independence-test-based d-separation determiner determines
%   d-separations by performing hypothesis tests of conditional
%   independence.
%
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.CIT.CITDSEPDETERMINER/ISDSEP
%   triggers a citPerformed event with
%   conditional-independence-test-performed data when a test is performed.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER,
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

% ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.CIT.CITDSEPDETERMINER
% implements the Adapter pattern.

properties(SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer

CITRCApplier % conditional-independence-test-reliability-criterion applier
CITPValueEstimator % conditional-independence-test-p-value estimator

type1ErrorRate % type I error rate -- a numeric real scalar in range [0, 1]
    
end

events

citPerformed % conditional-independence test performed

end

methods

% constructor

function Obj = citdsepdeterminer(CITRCApplier, CITPValueEstimator, alpha)
%CITDSEPDETERMINER Create conditional-independence-test-based d-separation determiner.
%   D = ORG.MENSXMACHINA.STATS.TESTS.CI.CITDSEPDETERMINER(A, E) creates a
%   conditional-independence-test-based d-separation determiner with
%	reliability-criterion applier A, p-value estimator E and type I error
%	rate 0.05. A and E have the same number of variables.
%
%   D = ORG.MENSXMACHINA.STATS.TESTS.CI.CITDSEPDETERMINER(A, E, ALPHA)
%   creates a determiner with type I error rate ALPHA. ALPHA is a numeric
%   real scalar in range [0, 1].

    % parse input
    assert(isa(CITRCApplier, 'citrcapplier') && isscalar(CITRCApplier));
    assert(isa(CITPValueEstimator, 'citpvalueestimator') && isscalar(CITPValueEstimator) && CITPValueEstimator.nVars == CITRCApplier.nVars);
    
    if nargin < 3
        alpha = 0.05;
    else
        validateattributes(alpha, {'numeric'}, {'real', 'scalar', '>=', 0, '<=', 1});
    end
    
    % set properties
    
    Obj.nVars = CITRCApplier.nVars;
    
    Obj.CITRCApplier = CITRCApplier;
    Obj.CITPValueEstimator = CITPValueEstimator;
    
    Obj.type1ErrorRate = alpha;

end

% abstract property implementations

function tf = isdsep(Obj, i, j, k)
    
    import citperformeddata;

    if isreliablecit(Obj.CITRCApplier, i, j, k) % test is reliable, attempt it

        % calculate test p-value
        [p stat] = citpvalue(Obj.CITPValueEstimator, i, j, k);

        % compare p-value with alpha
        tf = p > Obj.type1ErrorRate;

    else
        
        p = NaN;
        stat = NaN;

        tf = false;

    end
    
    notify(Obj, 'citPerformed', citperformeddata(i, j, k, p, stat));

end

function msc = maxsepsetcard(Obj, i, j, ind)
    
    % return CI test RC's best-case maximal conditioning set cardinality
    msc = bestmaxcondsetcard(Obj.CITRCApplier, i, j, ind);
    
end

end

end