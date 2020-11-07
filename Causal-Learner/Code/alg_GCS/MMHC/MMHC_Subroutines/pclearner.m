classdef pclearner < handle
%PCLEARNER Parents-and-children learner.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.PC.PCLEARNER is the abstract class of
%   parents-and-children learners. A parents-and-children learner learns,
%   for a set of random variables, the parents and children of a variable
%   in a Bayesian network representing the joint probability distribution
%   of the variables.

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

%LEARNPC Learn parents and children of a Bayesian network variable.
%   IND = LEARNPC(PCL, I), where PCL is a parents-and-children learner and
%   I is the linear index of PCL variable X, runs PCL with target X. I is a
%   numeric integer in range [1, M], where M is the number of PCL
%   variables. IND is a numeric row vector of distinct integers in range
%   [1, I - 1] or [I + 1, M]. Each element of IND is the linear index of a
%   presumed parent or child of X.
ind = learnpc(Obj, i);

end

methods(Access = protected)
    
% input parsers    

function parsecitpvalueinput(Obj, i)
%PARSELEARNPCINPUT Parse ORG.MENSXMACHINA.PGM.BN.LEARNING.PC.PCLEARNER/LEARNPC input.
%   PARSELEARNPCINPUT(PCL, ...), when PCL is a parents-and-children
%   learner, throws an error if its input is not valid input for
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.PC.PCLEARNER/LEARNPC.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.PC.PCLEARNER/LEARNPC.
    
    validateattributes(i, {'numeric'}, {'scalar', 'positive', '<=', Obj.nVars, 'integer'});
    
end

end

end