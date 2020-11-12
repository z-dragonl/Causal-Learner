classdef glglearner < lglearner
%GLGLEARNER Generalized-local-to-global learner.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.GLGLEARNER is the class of
%   generalized-local-to-global learners. Generalized-local-to-global
%   learners are both parent-and-children and skeleton learners.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.LGLEARNER,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.PC.PCLEARNER,
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

DSepDeterminer % d-separation determiner

maxSepsetCard % maximal sepset cardinality -- a numeric nonnegative integer
    
end

methods

% constructor

function Obj = glglearner(DSepDeterminer, varargin)
%GLGLEARNER Create generalized local learner.
%   L = ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.GLGLEARNER(D) creates
%   a generalized-local-to-global with d-separation determiner D.
%
%   [...] = ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.GLGLEARNER(D,
%   'Param1', VAL1, 'Param2', VAL2, ...) specifies additional parameter
%   name/value pairs chosen from the following:
%
%   'symCorrEnabled'    Whether symmetry correction is performed -- a
%                       logical scalar. Default is true.
%
%   'maxSepsetCard'     Maximal sepset cardinality -- a numeric integer in
%                       range [0 (M - 2)], where M is the number of
%                       variables of D.  Default is M - 2.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER.

    ip = inputParser;
    
    ip.addRequired('DSepDeterminer', @(a) isa(a, 'dsepdeterminer') && isscalar(a));
  
    ip.addParamValue('symCorrEnabled', true); % to be validated in LGLEARNER
    ip.addParamValue('maxSepsetCard', [], @(a) validateattributes(a, {'numeric'}, {'scalar', 'nonnegative', '<=', DSepDeterminer.nVars - 2, 'integer'}));
     
    ip.parse(DSepDeterminer, varargin{:});
    
    if isempty(ip.Results.maxSepsetCard)
        maxSepsetCard = DSepDeterminer.nVars - 2;
    else
        maxSepsetCard = ip.Results.maxSepsetCard;
    end
    
    % call ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.LGLEARNER
    % constructor
    Obj = Obj@lglearner(DSepDeterminer.nVars, ip.Results.symCorrEnabled);
    
    % set properties
    Obj.DSepDeterminer = DSepDeterminer;
    Obj.maxSepsetCard = maxSepsetCard;
    
end

end

methods(Sealed)

% abstract method implementations

function xTpcInd = learntpc(Obj, i)

    % create TPC(X), the set of tentative parents and children of X

    % initialize TPC(X), the set of tentative parents and children of X
    xTpcInd = initializetpc(Obj, i);
    
    % initialize OPEN(X) to V\({X} U TPC(X) U NTPC(X))

    % initialize OPEN(X) to V
    xOpenInd = 1:Obj.nVars;

    % remove X and TPC(X) from OPEN(X)
    xOpenInd = setdiff(xOpenInd, [i xTpcInd]);

    % init OPEN(X) priority queue
    xOpenPriorityQueue = inititializeopenpriorityqueue(Obj, xOpenInd);

    while true % apply interleaving strategy...
        
        % apply inclusion heuristic function

        % a) prioritize variables in OPEN(X) for inclusion in TPC(X) and
        % b) throw away non-eligible variables from OPEN(X)

%         fprintf('\nStep A and B...\n');
        
        [xOpenInd xOpenPriorityQueue] = ...
            updateopen(Obj, i, xTpcInd, xOpenInd, xOpenPriorityQueue);
                             
%         fprintf('\nOPEN(%d) = [%s]\n', i, num2str(xOpenInd));

        % c) insert in TPC(X) the highest priority variable(s) in OPEN(X) and
        % remove them from OPEN(X)
 
%         fprintf('\nStep C...\n');
        
        [xTpcInd xOpenInd xOpenPriorityQueue] = ...
            removefromopen(Obj, i, xTpcInd, xOpenInd, xOpenPriorityQueue);

%          fprintf('\nOPEN(%d) = [%s]\n', i, num2str(xOpenInd));
%          fprintf('\nTPC(%d) = [%s]\n', i, num2str(xTpcInd));
%          fprintf('\nStep D...\n');

        % apply elimination strategy to remove variables from TPC(X)
        xTpcInd = removefromtpc(Obj, i, xTpcInd, xOpenInd);
        
%         fprintf('\nTPC(%d) = [%s]\n', i, num2str(xTpcInd));

        if mustterminate(Obj, i, xTpcInd, xOpenInd)
            % ... until the termination criterion is met
            break;        
        end

    end

end

end

methods(Abstract, Access = protected)
    
% Note: Help for the following methods will be added in the future.

xTpcInd = initializetpc(Obj, xInd);  

xOpenPriorityQueue = inititializeopenpriorityqueue(Obj, xOpenInd);

[xOpenInd xOpenPriorityQueue] = ...
    updateopen(Obj, i, xTpdInd, xOpenInd, xOpenPriorityQueue);

[xTpcInd xOpenInd xOpenPriorityQueue] = ...
    removefromopen(Obj, i, xTpcInd, xOpenInd, xOpenPriorityQueue);

xTpcInd = removefromtpc(Obj, i, xTpcInd, xOpenInd);

tf = mustterminate(Obj, i, xTpcInd, xOpenInd);

end

end