classdef structurelearner < handle
%STRUCTURELEARNER Skeleton learner.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.STRUCTURELEARNER is the
%   abstract class of structure learners. A structure learner learns, for a
%   set of random variables, the structure of a Bayesian network
%   representing the joint probability distribution of the variables.

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

%LEARNSTRUCTURE Learn Bayesian network structure.
%   STRUCTURE = LEARNSTRUCTURE(SL, I) runs structure learner SL. SL is an
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.STRUCTURELEARNER Object.
%   STRUCTURE is a M-by-M sparse matrix. Each nonzero element in the lower
%   triangle of STRUCTURE denotes an edge in the structure.
structure = learnstructure(Obj);

end

end