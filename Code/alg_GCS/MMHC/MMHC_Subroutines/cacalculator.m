classdef cacalculator < handle
%CACALCULATOR Conditional-association calculator.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CACALCULATOR is the
%   abstract class of conditional-association calculator. A
%   conditional-association calculator calculates, for a set of variables,
%   a measure of association of a pair of variables given a subset of the
%   rest variables.

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

properties(Abstract, SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer
    
end

methods(Abstract)

%CA Calculate conditional association.
%   [A1 A2] = CA(C, I, J, K), where C is a conditional-association
%   calculator, I and J is the linear index of C variable X and Y,
%   respectivelly, K are the linear indices of subset Z of D variables, and
%   sets {X}, {Y} and Z are disjoint, calculates the value of a primary and
%   secondary association measure, A1 and A2, respectivelly, of X and Y
%   given Z. Secondary associations are used to break ties when comparing
%   associations.
%
%   Example:
%
%       % CACalculator is a conditional-association calculator
%
%       % calculate assoc(1, 2| [3 4])
%       [ca1 ca2] = ca(CACalculator, 1, 2, [3 4]);
[ca1 ca2] = ca(Obj, i, j, k);

end

methods(Access = protected)
    
% input parsers    

function parsecainput(Obj, i, j, k)
%PARSECAINPUT Parse ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CACALCULATOR/CA input.
%   PARSECAINPUT(D, ...), when D is a d-separation determiner, throws an
%   error if its input is not valid input for
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CACALCULATOR/CA.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CACALCULATOR/CA.
    
    validateattributes(i, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    
    validateattributes(j, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    
    validateattributes(k, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    
    u = [i j k];
    assert(length(unique(u)) == length(u));
    
end

end

end