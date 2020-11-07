classdef localscorer < handle
%LOCALSCORER Local scorer.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER is the
%   abstract class of local scorers. A local scorer scores, for a set of
%   random variables, the family of a variable (that is, the variable
%   itself and its parents) in a Bayesian network representing the joint
%   probability distribution of the variables.

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
%     Bayesian networks: The combination of knowledge and statistical data.
%     Machine Learning, 20, 197-243.

properties(Abstract, SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer
    
end

methods(Abstract)

%LOGLOCALSCORE Compute a family log score.
%   LS = LOGLOCALSCORE(LSR, I, J), where LSR is a local scorer, I is the
%   linear index of LSR variable X, and J are the linear indices of a set
%   PA(X) of LSR variables other than X, computes the log score of family
%   [X PA(X)]. LSR is a an
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER
%   Object. I is a numeric integer in range [1, M], where M is the number
%   of LSR variables. J is a numeric row vector of distinct integers in
%   range [1, I - 1] or [I + 1, M]. LS is numeric real scalar.
lsc = loglocalscore(Obj, i, j);

end

methods(Access = protected)
    
% input parsers    

function parsecitpvalueinput(Obj, i, j)
%PARSELOGLOCALSCOREINPUT Parse ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER/LOGLOCALSCORE input.
%   PARSELOGLOCALSCOREINPUT(LSR, ...), when LSR is an
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER array,
%   throws an error if its input is not valid input for
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER/LOGLOCALSCORE.
%
%   See also
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER/LOGLOCALSCORE.
    
    validateattributes(i, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    
    validateattributes(j, {'numeric'}, {'row', 'positive', '<=', Obj.nVars, 'integer'});
    assert(~ismember(i, j));
    
end

end

end