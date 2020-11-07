classdef lglearner < pclearner & skeletonlearner
%LGLEARNER Local-to-global learner.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.LGLEARNER is the class of
%   local-to-global learners. Local-to-global learners are both
%   parent-and-children and skeleton learners.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.PC.PCLEARNER,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.SKELETON.SKELETONLEARNER,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER.

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
% [1] Aliferis, C. F., Statnikov, A., Tsamardinos, I., Mani, S., and
%     Koutsoukos, X. D. 2010. Local Causal and Markov Blanket Induction for
%     Causal Discovery and Feature Selection for Classification Part I:
%     Algorithms and Empirical Evaluation. J. Mach. Learn. Res. 11 (Mar.
%     2010), 171-234.
% [2] I. Tsamardinos, L.E. Brown, C.F. Aliferis. "The Max-Min Hill-Climbing
%     Bayesian Network structure Learning Algorithm", Machine Learning
%     Journal; 65: 31-78

properties(SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer

symCorrEnabled % whether symmetry correction is performed -- a logical scalar
    
end

methods

% constructor

function Obj = lglearner(nVars, symCorrEnabled)
%LGLEARNER Create generalized local learner.
%   L = ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.LGLEARNER(N) creates a
%   local-to-global learner with N variables and symmetry correction
%   enabled.
%
%   L = ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.LGLEARNER(N, SC) also
%   specifies whether symmetry correction is enabled. SC is true when
%   symmetry correction is enabled and false otherwise.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER.

    validateattributes(nVars, {'numeric'}, {'scalar', 'nonnegative', 'integer'});
    
    if nargin < 2
        symCorrEnabled = true;
    else
        validateattributes(symCorrEnabled, {'logical'}, {'scalar'});
    end
    
    % set properties
    Obj.nVars = nVars;
    Obj.symCorrEnabled = symCorrEnabled;
end

end

methods(Sealed)

% abstract method implementations

function xPCInd = learnpc(Obj, i)

    import glglearner;

    % learn TPC(X)
%     fprintf('\nLearning TPC(%d)...\n', i);

    xTpcInd = learntpc(Obj, i);

    if Obj.symCorrEnabled
        
        % perform symmetry correction
        
        xPCInd = zeros(1, 0);

        for yInd = sort(xTpcInd) % for each variable Y in sorted TPC(X)

            % learn TPC(Y)
%             fprintf('\nLearning TPC(%d)...\n', yInd);

            yTpcInd = learntpc(Obj, yInd);
            
            if ismember(i, yTpcInd)
                xPCInd = [xPCInd yInd];
            end

        end
    
    else
        
        xPCInd = xTpcInd;
        
    end

end

function skeleton = learnskeleton(Obj)

    import glglearner;
    
    % initialize skeleton
    skeleton = sparse(Obj.nVars, Obj.nVars);

    for i = 1:Obj.nVars % for each variable
        
%         fprintf('\nLearning TPC(%d)...\n', i);

        % learn TPC(X)
        xTpcInd = learntpc(Obj, i);

        % update skeleton
        
        % update previously established links

        prevInd = 1:(i - 1);        
        prevTpc = xTpcInd(xTpcInd < i);        
        update = sparse(ones(size(prevTpc)), prevTpc, true(size(prevTpc)), 1, i - 1);
        
        if Obj.symCorrEnabled % symmetry correction on            
            skeleton(i, prevInd) = skeleton(i, prevInd) & update;
        else % symmetry correction off
            skeleton(i, prevInd) = skeleton(i, prevInd) | update;
        end   
        
        % establish new links
        skeleton(xTpcInd(xTpcInd > i), i) = true;

    end

end

end

methods(Abstract)

%LEARNTPC Learn tentative parents and children of a Bayesian network variable.
%   IND = LEARNTPC(LGL, I), where LGL is a local-to-global learner and I is
%   the linear index of LGL variable X, learns a set of tentative parents
%   and children of X, a superset of the set of parents and children of X.
%   I is a numeric integer in range [1, M], where M is the number of LGL
%   variables. IND is a numeric row vector of distinct integers in range
%   [1, I - 1] or [I + 1, M]. Each element of IND is the linear index of a
%   tentative parent or child of X.
xTpcInd = learntpc(Obj, i)

end

end