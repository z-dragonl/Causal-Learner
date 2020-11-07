classdef citrcapplier < handle
%CITRCAPPLIER Conditional-independence-test-reliability-criterion applier.
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER is the class of
%   conditional-independence-test-reliability-criterion appliers. A
%   conditional-independence-test-reliability-criterion applier determines,
%   for a set of variables, whether hypothesis tests of conditional
%   independence involving variables of the set are reliable according to a
%   reliability criterion. A
%   conditional-independence-test-reliability-criterion applier also bounds
%   the maximal conditioning-set cardinality for which tests are reliable.

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

%ISRELIABLECIT Determine whether a hypothesis test of conditional independence is reliable.
%   TF = ISRELIABLECIT(CITRCAPPOBJ, I, J, K), where CITRCAPPOBJ is a
%   conditional-independence-test-reliability-criterion applier, I and J is
%   the linear index of CITRCAPPOBJ variable X and Y, respectivelly, and K
%   are the linear indices of set Z of CITRCAPPOBJ variables, returns
%   logical 1 (true) if the test of the hypothesis that X and Y are
%   independent given S is reliable according to the reliavility criterion
%   of CITRCAPPOBJ and logical 0 (false) otherwise. I and J are numeric
%   integers in range [1, M], where M is the number of CITRCAPPOBJ
%   variables. K is a numeric row vector of integers in range [1, M].
%
%   Example:
%
%       % rcApplier is a %
%       conditional-independence-test-reliability-criterion applier
%
%       % determine whether the testing of Ind(1, 2| [3 4]) is reliable tf
%       = isreliablecit(rcApplier, 1, 2, [3 4]);
tf = isreliablecit(Obj, i, j, k);

%WORSTMAXCONDSETCARD Lower-bound worst-case maximal conditioning-set cardinality.
%   LB = WORSTMAXCONDSETCARD(CITRCAPPLIER, I, J, IND), where CITRCAPPLIER
%   is a conditional-independence-test-reliability-criterion applier, I and
%   J is the linear index of CITRCAPPLIER variable X and Y, respectivelly,
%   and IND are the linear indices of set S of CITRCAPPLIER variables,
%   returns a lower bound on the XY worst-case maximal conditioning-set
%   cardinality (XY-worst-max-k) of S. The XY-worst-max-k of S is the
%   maximal K such that all tests of conditional independence of X and Y
%   given a subset of S with cardinality K are reliable according to the
%   reliability criterion of CITRCAPPLIER. I and J are numeric integers in
%   range [1, M], where M is the number of CITRCAPPOBJ variables. IND is a
%   numeric row vector of integers in range [1, M]. LB is a numeric integer
%   in range [0, LENGTH(IND)].
%
%   LB = WORSTMAXCONDSETCARD(CITRCAPPLIER, I, ZEROS(1, 0), IND) returns a
%   lower bound on the X worst-case maximal conditioning-set cardinality
%   (X-worst-max-k) of S. The X-worst-max-k of S is the maximal K such that
%   all tests of conditional independence of X and any variable Y of S
%   given a subset of S\{Y} with cardinality K are reliable according to
%   the reliability criterion of CITRCAPPLIER. LB is in range [0,
%   LENGTH(IND) - 1].
%
%   LB = WORSTMAXCONDSETCARD(CITRCAPPLIER, ZEROS(1, 0), I, IND) is the same
%   as LB = WORSTMAXCONDSETCARD(CITRCAPPLIER, I, ZEROS(1, 0), IND).
%
%   LB = WORSTMAXCONDSETCARD(CITRCAPPLIER, ZEROS(1, 0), ZEROS(1, 0), IND)
%   returns a lower bound on the worst-case maximal conditioning-set
%   cardinality (worst-max-k) of S. The worst-max-k of S is the maximal K
%   such that all tests of conditional independence of variables X and Y of
%   S given a subset of S\{X, Y} with cardinality K are reliable according
%   to the reliability criterion of CITRCAPPLIER. LB is in range [0,
%   LENGTH(IND) - 2].
lb = worstmaxcondsetcard(Obj, i, j, ind);

%BESTMAXCONDSETCARD Upper-bound best-case maximal conditioning-set cardinality.
%   UB = BESTMAXCONDSETCARD(CITRCAPPLIER, I, J, IND), where CITRCAPPLIER is
%   a conditional-independence-test-reliability-criterion applier, I and J
%   is the linear index of CITRCAPPLIER variable X and Y, respectivelly,
%   and IND are the linear indices of set S of CITRCAPPLIER variables,
%   returns a lower bound on the XY best-case maximal conditioning-set
%   cardinality (XY-best-max-k) of S. The XY-best-max-k of S is the maximal
%   K such that at least one test of conditional independence of X and Y
%   given a subset of S with cardinality K is reliable according to the
%   reliability criterion of CITRCAPPLIER. I and J are numeric integers in
%   range [1, M], where M is the number of CITRCAPPOBJ variables. IND is a
%   numeric row vector of integers in range [1, M]. UB is a numeric integer
%   in range [0, LENGTH(IND)].
%
%   UB = BESTMAXCONDSETCARD(CITRCAPPLIER, I, ZEROS(1, 0), IND) returns an
%   upper bound on the X best-case maximal conditioning-set cardinality
%   (X-best-max-k) of S. The X-best-max-k of S is the maximal K such that
%   at least one test of conditional independence of X and any variable Y
%   of S given a subset of S\{Y} with cardinality K is reliable according
%   to the reliability criterion of CITRCAPPLIER. UB is in range [0,
%   LENGTH(IND) - 1].
%
%   UB = BESTMAXCONDSETCARD(CITRCAPPLIER, ZEROS(1, 0), I, IND) is the same
%   as UB = BESTMAXCONDSETCARD(CITRCAPPLIER, I, ZEROS(1, 0), IND).
%
%   UB = BESTMAXCONDSETCARD(CITRCAPPLIER, ZEROS(1, 0), ZEROS(1, 0), S)
%   returns an upper bound on the best-case maximal conditioning-set
%   cardinality (best-max-k) of S. The best-max-k of S is the maximal K
%   such that at least one test of conditional independence of variables X
%   and Y of S given a subset of S\{X, Y} with cardinality K is reliable
%   according to the reliability criterion of CITRCAPPLIER. UB is in range
%   [0, LENGTH(IND) - 2].
ub = bestmaxcondsetcard(Obj, i, j, ind);

end

methods(Access = protected)
    
% input parsers    

function parseisreliablecitinput(Obj, i, j, k)
%PARSEISRELIABLECITINPUT Parse ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/ISRELIABLECIT input.
%   PARSEISRELIABLECITINPUT(CITRCAPPOBJ, ...), when CITRCAPPOBJ is a
%   conditional-independence-test-reliability-criterion applier, throws an
%   error if its input is not valid input for
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/ISRELIABLECIT.
%
%   See also ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/ISRELIABLECIT.
    
    validateattributes(i, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    validateattributes(j, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    validateattributes(k, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    
end

function parsemaxcondsetinput(Obj, i, j, ind)
%PARSEMAXCONDSETCARDINPUT Parse ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/WORSTMAXCONDSETCARD or ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/BESTMAXCONDSETCARD input.
%   PARSEMAXCONDSETCARDINPUT(CITRCAPPOBJ, ...), when CITRCAPPOBJ is a
%   conditional-independence-test-reliability-criterion applier, throws an
%   error if its input is not valid input for
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/WORSTMAXCONDSETCARD or
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/BESTMAXCONDSETCARD.
%
%   See also
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/WORSTMAXCONDSETCARD,
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CITRCAPPLIER/BESTMAXCONDSETCARD.
    
    validateattributes(i, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    assert(length(i) < 1);
    
    validateattributes(j, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    assert(length(j) < 1);
    
    validateattributes(ind, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});

end

end

end