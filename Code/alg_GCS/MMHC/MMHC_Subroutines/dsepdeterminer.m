classdef dsepdeterminer < handle
%DSEPDETERMINER D-separation determiner.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER is the abstract
%   class of d-separation determiners. A d-separation determiner determines
%   d-separation relationships among a set of variables.

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

properties(Abstract, SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer
    
end

methods(Abstract)

%ISDSEP Determine d-separation.
%   TF = ISDSEP(D, I, J, K), where D is a d-separation determiner, I and J
%   is the linear index of D variable X and Y, respectivelly, K are the
%   linear indices of subset Z of D variables, and sets {X}, {Y} and Z are
%   disjoint, returns true if X and Y are d-separated by Z and false
%   otherwise.
%
%   Example:
%
%       % DSepDeterminer is a d-separation determiner
%
%       % determine if DSep(1, 2| [3 4])
%       tf = isdsep(DSepDeterminer, 1, 2, [3 4]);
tf = isdsep(Obj, i, j, k);

%MAXSEPSETCARD Upper-bound sepset cardinality.
%   C = MAXSEPSETCARD(D, I, J, IND), where D is a d-separation determiner,
%   I and J is the linear index of D variable X and Y, respectivelly, IND
%   are the linear indices of set S of D variables, and sets {X}, {Y} and S
%   are disjoint, returns a lower bound on the XY maximal
%   sepset-cardinality (XY-max-sc) of S. The XY-max-sc of S is the maximal
%   cardinality of a subset of S that could d-separate X and Y. I and J are
%   numeric integers in range [1, M], where M is the number of D variables.
%   IND is a numeric row vector of integers in range [1, M]. C is a numeric
%   integer in range [0, LENGTH(IND)].
%
%   C = MAXSEPSETCARD(D, I, ZEROS(1, 0), IND) returns an upper bound on the
%   X maximal sepset-cardinality (X-max-sc) of S. The X-max-sc of S is the
%   maximal cardinality of a subset of S\{Y} that could d-separate X from
%   some variable Y of S. C is in range [0, LENGTH(IND) - 1].
%
%   C = MAXSEPSETCARD(D, ZEROS(1, 0), I, IND) is the same as C =
%   MAXSEPSETCARD(D, I, ZEROS(1, 0), IND).
%
%   C = MAXSEPSETCARD(D, ZEROS(1, 0), ZEROS(1, 0), S) returns an upper
%   bound on the maximal sepset-cardinality (max-sc) of S. The max-sc of S
%   is the maximal cardinality of a subset of S\{X, Y} that could
%   d-separate distinct variables X and Y of S. C is in range [0,
%   LENGTH(IND) - 2].
msc = maxsepsetcard(Obj, i, j, ind);

end

methods(Access = protected)
    
% input parsers    

function parseisdsepinput(Obj, i, j, k)
%PARSEISDSEPINPUT Parse ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER/ISDSEP input.
%   PARSEISDSEPINPUT(D, ...), when D is a d-separation determiner, throws
%   an error if its input is not valid input for
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER/ISDSEP.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER/ISDSEP.
    
    validateattributes(i, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    
    validateattributes(j, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    
    validateattributes(k, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    
    u = [i j k];
    assert(length(unique(u)) == length(u));
    
end

function parsemaxsepsetcardinput(Obj, i, j, ind)
%PARSEMAXSEPSETCARDINPUT Parse ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER/MAXSEPSETCARD input.
%   PARSEMAXSEPSETCARDINPUT(D, ...), when D is d-separation determiner,
%   throws an error if its input is not valid input for
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER/MAXSEPSETCARD.
%
%   See also
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER/MAXSEPSETCARD.
    
    validateattributes(i, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    assert(length(i) < 1);
    
    validateattributes(j, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    assert(length(j) < 1);
    
    validateattributes(ind, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    
    u = [i j ind];
    assert(length(unique(u)) == length(u));

end

end

end