classdef(Sealed) hillclimber < structurelearner
%HILLCLIMBER Hill climber.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.HS.HILLCLIMBER is
%   the class of hill climbers.

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
    
end

properties(GetAccess = private, SetAccess = immutable)

LocalScorer % local scorer -- an ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER Object

maxNBarrenRepeats % maximal number of iterations without increase in maximal score -- a numeric nonnegative integer
tabuListLength % TABU-list length -- a numeric nonnegative integer

initialStructure % initial structure -- an M-by-M sparse matrix
candidateParentMatrix % candidate-parent matrix -- either an M-by-M sparse matrix or []

end

methods

function Obj = hillclimber(LocalScorer, varargin)
%HILLCLIMBER Create hill climber.
%   H =
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.HS.HILLCLIMBER(S)
%   creates a hill climber with local scorer S. S is an
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER
%   Object.
%
%   H =
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.HS.HILLCLIMBER(S,
%   'Param1', VAL1, 'Param2', VAL2, ...) specifies additional parameter
%   name/value pairs chosen from the following:
%
%       'candidateParentMatrix'     Candidate-parent matrix -- either a
%                                   sparse M-by-M sparse matrix, where M is
%                                   the number of variables of S, or []
%                                   (default). In the first case, a nonzero
%                                   or true CANDIDATEPARENTMATRIX(bestActionXInd, bestActionYInd)
%                                   indicates that node bestActionXInd is a candidate
%                                   parent of node bestActionYInd. In the second case,
%                                   all other nodes are considered to be
%                                   candidate parents of a node.
%
%       'initialStructure'          Initial structure -- a sparse M-by-M
%                                   sparse matrix, where M is the number of
%                                   variables of S. Each nonzero element in
%                                   the lower triangle of the matrix
%                                   denotes an edge in the structure.
%                                   Default is the M-by-M all-zero sparse
%                                   matrix.
%
%       'maxNBarrenRepeats'         Maximal number of iterations
%                                   without increase in maximal score
%                                   before the algorithm stops -- a numeric
%                                   nonnegative integer. Default is 15.
%
%       'tabuListLength'            TABU list length -- a numeric
%                                   nonnegative integer. Default is 100.
%
%   See also
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.STRUCTURE.SNS.LOCAL.LOCALSCORER.




    ip = inputParser;
    
    ip.addRequired('LocalScorer', @(a) isa(a, 'localscorer'));
    
    ip.addParamValue('maxNBarrenRepeats', 15, @(a) validateattributes(a, {'numeric'}, {'nonnegative', 'integer', 'scalar'}));
    ip.addParamValue('tabuListLength', 100, @(a) validateattributes(a, {'numeric'}, {'nonnegative', 'integer', 'scalar'}));
    ip.addParamValue('initialStructure', [], @(a) issparse(a) && size(a, 1) == LocalScorer.nVars && size(a, 2) == LocalScorer.nVars);
    ip.addParamValue('candidateParentMatrix', [], @(a) issparse(a) && size(a, 1) == LocalScorer.nVars && size(a, 2) == LocalScorer.nVars);


    ip.parse(LocalScorer, varargin{:});

    nVars = LocalScorer.nVars;
    
    if isempty(ip.Results.initialStructure)
        initialStructure = logical(sparse(nVars, nVars));
    else
        initialStructure = logical(ip.Results.initialStructure);
    end

    if isempty(ip.Results.candidateParentMatrix)
        candidateParentMatrix = [];
    else
        candidateParentMatrix = logical(ip.Results.candidateParentMatrix);
    end
    
    % set properties
    Obj.nVars = LocalScorer.nVars;
    Obj.LocalScorer = LocalScorer;
    Obj.maxNBarrenRepeats = ip.Results.maxNBarrenRepeats;
    Obj.tabuListLength = ip.Results.tabuListLength;
    Obj.initialStructure = initialStructure;
    Obj.candidateParentMatrix = candidateParentMatrix;
    
end

% abstract method implementations

function [structure logScore localLogScores] = learnstructure(Obj)
%LEARNSTRUCTURE Learn Bayesian network structure.
%   STRUCTURE = LEARNSTRUCTURE(HC, I) runs hill climber HC. STRUCTURE is a
%   M-by-M sparse matrix. Each nonzero element in the lower triangle of
%   STRUCTURE denotes an edge in the structure.
%
%   [STRUCTURE LOGSCORE] = LEARNSTRUCTURE(HC, I) also returns the log score
%   of STRUCTURE. LOGSCORE is a numeric real scalar.
%
%   [STRUCTURE LOGSCORE LOCALLOGSCORES] = LEARNSTRUCTURE(HC, I) also
%   returns the log score of each node family in STRUCTURE. LOCALLOGSCORES
%   is an 1-by-M numeric real scalars. Each element of LOCALLOGSCORES is
%   the log score of the corresponding family in STRUCTURE. The elements of
%   LOCALLOGSCORES sum to LOGSCORE.
    
%     import org.mensxmachina.pgm.bn.learning.structure.sns.local.hc.hillclimber;

    % determine whether the hill climber is constrained
    constrained = ~isempty(Obj.candidateParentMatrix);
    

    if constrained % constrained HC

        % convert candidate-parent matrix to vector
        cpmVec = Obj.candidateParentMatrix(:);

        % find candidate parent indices in CPM vector
        cpIndInCpmVec = find(cpmVec);

    else % plain HC

        vInd = 1:Obj.nVars; % create the set of nodes

    end
    

    % initialize structure
    structure = Obj.initialStructure;

    % get structure size
    structureSize = size(structure);
    

    % initialize local log scores
    
    localLogScores = zeros(1, Obj.nVars);
    
    for iVar = 1:Obj.nVars
        localLogScores(iVar) = Obj.LocalScorer.loglocalscore(iVar, find(structure(:, iVar))');
    end
    

    % sum local log scores to obtain structure log-score
    logScore = sum(localLogScores);
    

    % initialize ancestor matrix
    ancestors = ancestormatrix(structure);

    % initialize successor log-score differences and local log scores

    deltaLogScore = cell(1, 3);
    
    xFAScoreInSuccessor = cell(1, 3);
    yFAScoreInSuccessor = cell(1, 3);

    nzmax = nnz(Obj.candidateParentMatrix);
    
    m = Obj.nVars*Obj.nVars;

    for iAction = 1:3 % for each action
        
        deltaLogScore{iAction} = spalloc(m, 1, nzmax);
        
        xFAScoreInSuccessor{iAction} = spalloc(m, 1, nzmax);
        yFAScoreInSuccessor{iAction} = spalloc(m, 1, nzmax);
        
    end

    % update log-score differences
    
    for iVar = 1:Obj.nVars % for each node
        localupdatedeltalogscores(iVar);
    end

    % initialize TABU list
    
    if Obj.tabuListLength > 0 % a TABU list is used
        
        tabu = cell(1, Obj.tabuListLength); 
        tabu{1} = structure;
        
    end

    % initialize counters
    iRepeat = 0; % iteration counter
    nBarrenRepeats = 0; % barren-iteration counter

    % initialize best structure
    bestStructure = structure;
    bestStructureLogScore = logScore;
    bestStructureLocalLogScores = localLogScores;

    while true

        % update repetition counter
        iRepeat = iRepeat + 1;
        

        % get best successor not in TABU list

        % copy log-score differences
        workingDeltaLogScore = deltaLogScore;

        aSuccessorIsFound = false;

%         if ip.Results.Verbose
%             fprintf('\n\nSelecting successor...\n');
%         end

        while ~aSuccessorIsFound % while a successor is not found

            actionMaxDeltaLogScoreXInd = zeros(1, 3);
            actionMaxDeltaLogScoreYInd = zeros(1, 3);
            actionMaxDeltaLogScoreInd = zeros(1, 3);
            actionMaxDeltaLogScore = zeros(1, 3);

            for iAction = 1:3 % for each bestActionInd

                % find best parameters

                if constrained % constrained HC

                    % find best parameters' maxDeltaLogScore and index among candidate parents 
                    [actionMaxDeltaLogScore(iAction) thisActionMaxDeltaLogScoreInd] = ...
                        max(workingDeltaLogScore{iAction}(cpIndInCpmVec));

                    % find index among all possible parents
                    actionMaxDeltaLogScoreInd(iAction) = ...
                        cpIndInCpmVec(thisActionMaxDeltaLogScoreInd);   

                else % HC

                    % find best parameters' maxDeltaLogScore and index
                    [actionMaxDeltaLogScore(iAction) actionMaxDeltaLogScoreInd(iAction)] = ...
                        max(workingDeltaLogScore{iAction}); 

                end

                % convert to subscripts
                [actionMaxDeltaLogScoreXInd(iAction) actionMaxDeltaLogScoreYInd(iAction)] = ...
                    ind2sub(structureSize, actionMaxDeltaLogScoreInd(iAction));

            end

            % find best action
            [maxDeltaLogScore bestActionInd] = max(actionMaxDeltaLogScore);

            if isnan(maxDeltaLogScore) % no such action    
                
%                 if ip.Results.Verbose
%                     fprintf('\n\nNo best action.');
%                 end

                break;
                
            end

            % get best bestActionInd parameters
            bestActionXInd = actionMaxDeltaLogScoreXInd(bestActionInd);
            bestActionYInd = actionMaxDeltaLogScoreYInd(bestActionInd);
            bestActionXYInd = actionMaxDeltaLogScoreInd(bestActionInd);
    
%             if ip.Results.Verbose
%                 fprintf('\nSelected ');
%                 if bestActionInd == 1
%                     fprintf('ADD(%d->%d)', bestActionXInd, bestActionYInd);
%                 elseif bestActionInd == 2
%                     fprintf('DEL(%d->%d)', bestActionXInd, bestActionYInd); 
%                 else
%                     fprintf('REV(%d->%d)', bestActionXInd, bestActionYInd);
%                 end
%                 fprintf(' (maxDeltaLogScore = %.5f)', maxDeltaLogScore);
%             end

            % apply best action
            
            if bestActionInd == 1
                
                % add X -> Y
                structure(bestActionXInd, bestActionYInd) = 1;
                
            elseif bestActionInd == 2
                
                % delete X -> Y
                structure(bestActionXInd, bestActionYInd) = 0;
                
            else
                
                % reverse X -> Y
                structure(bestActionXInd, bestActionYInd) = 0;
                structure(bestActionYInd, bestActionXInd) = 1;
                
            end

            if successorisdag() % no cycles
                
                aSuccessorIsFound = true;
                
                for iTabuStructure = 1:length(tabu) % for each structure in TABU list
                    
                    if isequal(structure, tabu{iTabuStructure}) 
                        
%                         if ip.Results.Verbose
%                             fprintf(' In TABU list.');
%                         end

                        aSuccessorIsFound = false;
                        
                        break;
                        
                    end
                    
                end
                
%             elseif ip.Results.Verbose
%                 fprintf(' Introduces cycles.');

            end

            if ~aSuccessorIsFound

                % reverse bestActionInd
                
                if bestActionInd == 1 % added X -> Y
                    
                    % delete X -> Y
                    structure(bestActionXInd, bestActionYInd) = 0; 
                    
                elseif bestActionInd == 2 % deleted X -> Y
                    
                    % add X -> Y
                    structure(bestActionXInd, bestActionYInd) = 1; 
                    
                else % reversed X -> Y
                    
                    % reverse Y -> X
                    structure(bestActionYInd, bestActionXInd) = 0;
                    structure(bestActionXInd, bestActionYInd) = 1;
                    
                end

                % set bestActionInd as unavailable
                workingDeltaLogScore{bestActionInd}(bestActionXYInd) = NaN;

            end

        end

        if ~aSuccessorIsFound % no such successor
            
%             fprintf('\n\nNo successor. Stopping...\n');
            
            break;
            
        end

        % update local log scores
        
        localLogScores(bestActionYInd) = yFAScoreInSuccessor{bestActionInd}(bestActionXInd, bestActionYInd);
        
        if bestActionInd == 3
            localLogScores(bestActionXInd) = xFAScoreInSuccessor{bestActionInd}(bestActionXInd, bestActionYInd);
        end
        
        % update structure log-score
        logScore = sum(localLogScores);

        if Obj.tabuListLength > 0 % a TABU list is used
            
            % insert structure in TABU list
            tabu{mod(iRepeat, Obj.tabuListLength) + 1} = structure; 
            
        end

%         if ip.Results.Verbose
%             fprintf('\n');
%         end

        % print bestActionInd
        
%         fprintf('\n%d) ', iRepeat);
%         
%         if bestActionInd == 1
%             fprintf('ADD(%d -> %d)', bestActionXInd, bestActionYInd);
%         elseif bestActionInd == 2
%             fprintf('DEL(%d -> %d)', bestActionXInd, bestActionYInd); 
%         else
%             fprintf('REV(%d -> %d)', bestActionXInd, bestActionYInd);
%         end
%         
%         fprintf(' (log score = %.5f', logScore);

        if logScore > bestStructureLogScore && abs(logScore - bestStructureLogScore) > 1e-10 % current structure scores higher than the current best structure

            % set current structure as the current best structure
            bestStructure = structure;
            bestStructureLogScore = logScore; 
            bestStructureLocalLogScores = localLogScores;

%             fprintf(', max)');

            nBarrenRepeats = 0; % reset stop counter

        else

            nBarrenRepeats = nBarrenRepeats + 1; % increase stop counter

%             fprintf(')');

        end

        if nBarrenRepeats == Obj.maxNBarrenRepeats % stop counter reached upper limit
            
%             fprintf('\n\n%d repetitions occured without increase in maximal score. Stopping...', nBarrenRepeats);
            
            break;
            
        end

        % update ancestors graph
        updateancestors();

%         assert(isequal(ancestors, org.mensxmachina.graph.ancestormatrix(structure)));

        % update successor scores
        updatedeltalogscores();

    end

    % set best structure as the output structure
    structure = bestStructure;
    logScore = bestStructureLogScore;
    localLogScores = bestStructureLocalLogScores;

    % nested functions

    function updatedeltalogscores()

        if bestActionInd < 3 % added or deleted X -> Y
            localupdatedeltalogscores(bestActionYInd);
        else % reversed X -> Y
            localupdatedeltalogscores(bestActionYInd);
            localupdatedeltalogscores(bestActionXInd);
        end

    end

    function localupdatedeltalogscores(i)

        if ~constrained % plain Hill Climbing

            for j = vInd % for each node

                localupdatecpdeltalogscores(i, j);
                localupdateccdeltalogscores(i, j);

            end

        else % constrained HC

%             if ip.Results.Verbose
%                 fprintf('\n\nEvaluating ACT(j->%d)...\n', i);
%             end

            for j = find(Obj.candidateParentMatrix(:, i))' % for each candidate parent of X
                localupdatecpdeltalogscores(i, j);
            end

%             if ip.Results.Verbose
%                 fprintf('\n\nEvaluating REV(%d->j)...\n', i);
%             end

            for j = find(Obj.candidateParentMatrix(i, :)) % for each candidate child of X
                localupdateccdeltalogscores(i, j);
            end

        end

    end

    function localupdatecpdeltalogscores(i, j)

        % get index of Y->X
        yxInd = sub2ind(structureSize, j, i);

%         if ip.Results.Verbose
%             fprintf('\nEvaluating ADD(%d->%d)... ', j, i);
%         end

        if j == i || structure(j, i) || structure(i, j) % either Y == X or Y->X or X->Y exists

            deltaLogScore{1}(yxInd) = NaN;

        else % neither Y == X nor Y->X or X->Y exists

            % add Y->X
            structure(j, i) = 1;

            [xFAScoreInSuccessor{1}(j, i) yFAScoreInSuccessor{1}(j, i)] = loglocalscoreinsuccessor(1, j, i);

            % attention: order of operations MATTERS!
            deltaLogScore{1}(yxInd) = xFAScoreInSuccessor{1}(j, i) - localLogScores(j) + yFAScoreInSuccessor{1}(j, i) - localLogScores(i);

            % delete Y->X
            structure(j, i) = 0;

        end

%         if ip.Results.Verbose
%             fprintf('Delta = %.5f ', deltaLogScore{1}(yxInd));
%             fprintf('\nEvaluating DEL(%d->%d)... ', j, i);
%         end

        if j == i || ~structure(j, i) % either Y == X or Y->X does not exist

            deltaLogScore{2}(yxInd) = NaN;

        else % Y->X exists

            structure(j, i) = 0; % delete Y->X

            [xFAScoreInSuccessor{2}(j, i) yFAScoreInSuccessor{2}(j, i)] = loglocalscoreinsuccessor(2, j, i);

            deltaLogScore{2}(yxInd) = xFAScoreInSuccessor{2}(j, i) - localLogScores(j) + yFAScoreInSuccessor{2}(j, i) - localLogScores(i);

            structure(j, i) = 1; % add Y->X

        end

%         if ip.Results.Verbose
%             fprintf('Delta = %.5f ', deltaLogScore{2}(yxInd));
%             fprintf('\nEvaluating REV(%d->%d)... ', j, i);
%         end

        % evaluate REV(Y->X)

        if j == i || ~structure(j, i) % either Y == X or Y->X does not exist

            deltaLogScore{3}(yxInd) = NaN;

        else % Y->X exists

            structure(j, i) = 0; % delete Y->X
            structure(i, j) = 1; % add X->Y

            [xFAScoreInSuccessor{3}(j, i) yFAScoreInSuccessor{3}(j, i)] = loglocalscoreinsuccessor(3, j, i);

            deltaLogScore{3}(yxInd) = xFAScoreInSuccessor{3}(j, i) - localLogScores(j) + yFAScoreInSuccessor{3}(j, i) - localLogScores(i);

            structure(j, i) = 1; % add Y->X
            structure(i, j) = 0; % delete X->Y

        end

%         if ip.Results.Verbose
%             fprintf('Delta = %.5f ', deltaLogScore{3}(yxInd));
%         end

    end

    function localupdateccdeltalogscores(i, j)

        % get index of X->Y
        xyInd = sub2ind(structureSize, i, j);    

%         if ip.Results.Verbose
%             fprintf('\nEvaluating REV(%d->%d)... ', i, j);
%         end

        if i == j || ~structure(i, j) % either X == Y or X->Y does not exist

            deltaLogScore{3}(xyInd) = NaN;

        else % X->Y exists

            structure(i, j) = 0; % delete X->Y
            structure(j, i) = 1; % add Y->X

            [xFAScoreInSuccessor{3}(i, j) yFAScoreInSuccessor{3}(i, j)] = loglocalscoreinsuccessor(3, i, j);

            deltaLogScore{3}(xyInd) = yFAScoreInSuccessor{3}(i, j) - localLogScores(i) + xFAScoreInSuccessor{3}(i, j) - localLogScores(j);

            structure(i, j) = 1; % add X->Y
            structure(j, i) = 0; % delete Y->X

        end

%         if ip.Results.Verbose
%             fprintf('Delta = %.5f ', deltaLogScore{3}(xyInd));
%         end

    end

    function updateancestors()

        if bestActionInd == 1 % added X -> Y
            
            updateancestorsonadd(bestActionXInd, bestActionYInd);
            
        elseif bestActionInd == 2 % deleted X -> Y
            
            updateancestorsondelete(bestActionYInd); 
            
        else % reversed X -> Y
            
            % delete Y -> X
            structure(bestActionYInd, bestActionXInd) = 0;
            
            updateancestorsondelete(bestActionYInd);
            
            % add Y -> X
            structure(bestActionYInd, bestActionXInd) = 1;
            
            updateancestorsonadd(bestActionYInd, bestActionXInd);
            
        end

    end

    function updateancestorsonadd(i, j)

        % set X as an ancestor of Y
        ancestors(i, j) = 1;

        % set ancestors of X as ancestors of Y
        ancestors(full(ancestors(:, i)), j) = 1;

        % set ancestors of X as ancestors of the descendants of X

        for iXDE = find(ancestors(j, :)) % for each descendant of X

            % set ancestors of X as ancestors of this descendant of X
            ancestors(full(ancestors(:, j)), iXDE) = 1;

        end

    end

    function updateancestorsondelete(j)

        % set ancestors of Y to parents of Y
        ancestors(:, j) = structure(:, j);

        % set ancestors of parents of Y to ancestors of Y

        for iYPA = find(structure(:, j))' % for each parent of Y

            % set ancestors of this parent as ancestors of Y
            ancestors(full(ancestors(:, iYPA)), j) = 1;

        end

        % find descendants of Y
        yDE = find(ancestors(j, :));

        % get subgraph of descendants of Y
        yDEStructure = structure(yDE, yDE);

        if ~issparse(yDEStructure)
            yDEStructure = sparse(yDEStructure);
        end

        % find topological order
        order = graphtopoorder(yDEStructure);

        for iYDE = yDE(order) % for each descendant of Y in topological order

            % set ancestors of this descendant to parents of this descendant
            ancestors(:, iYDE) = structure(:, iYDE);

            % set ancestors of parents of this descendant to ancestors of
            % this descendant

            for thisXDEIPA = find(structure(:, iYDE))' % for each parent of this descendant

                % set ancestors of this parent of this descendant as
                % ancestors of this descendant
                ancestors(full(ancestors(:, thisXDEIPA)), iYDE) = 1;

            end

        end

    end
    
    function tf = successorisdag()

        if bestActionInd == 1 % added X -> Y

            % DAG iff Y is not an ancestor of Y in current structure
            tf = ~ancestors(bestActionYInd, bestActionXInd);

        elseif bestActionInd == 2 % deleted X -> Y

            tf = true;

        else % 3, reversed X -> Y

            % structure iff X is not an ancestor (in current structure) of
            % a parent (in successor structure) of Y

            % (X is not its ancestor anyway so we could also use the
            % parents in the current structure)

            tf = ~any(ancestors(bestActionXInd, full(structure(:, bestActionYInd))));

        end

    end

    function [logXFaScoreInSuccessor logYFaScoreInSuccessor] = loglocalscoreinsuccessor(actionInd, actionXInd, actionYInd)

        if actionInd == 1 % add X -> Y

            logXFaScoreInSuccessor = localLogScores(actionXInd);
            logYFaScoreInSuccessor = Obj.LocalScorer.loglocalscore(actionYInd, find(structure(:, actionYInd))');

        elseif actionInd == 2 % 2, delete X -> Y

            logXFaScoreInSuccessor = localLogScores(actionXInd);
            logYFaScoreInSuccessor = Obj.LocalScorer.loglocalscore(actionYInd, find(structure(:, actionYInd))');

        else % reverse X -> Y

            logXFaScoreInSuccessor = Obj.LocalScorer.loglocalscore(actionXInd, find(structure(:, actionXInd))');
            logYFaScoreInSuccessor = Obj.LocalScorer.loglocalscore(actionYInd, find(structure(:, actionYInd))');

        end

    end

end

end

end