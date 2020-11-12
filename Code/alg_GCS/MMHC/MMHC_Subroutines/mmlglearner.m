classdef(Sealed) mmlglearner < glglearner
%MMLGLEARNER Max-min local-to-global learner.
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.MMLGLEARNER is the class
%   of max-min local-to-global learners. Max-min local-to-global learners
%   are both parent-and-children and skeleton learners.
%
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.MMLGLEARNER/LEARNPC and
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.MMLGLEARNER/LEARNSKELETON
%   trigger a sepsetIdentified event with sepset-identified data when a
%   sepset is identified.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.GLGLEARNER,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.PC.PCLEARNER,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.SKELETON.SKELETONLEARNER,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CACALCULATOR,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.SEPSETIDENTIFIEDDATA.

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
    
CACalculator % conditional-association calculator
    
end

events
    
sepsetIdentified % sepset identified

end

methods

% constructor

function Obj = mmlglearner(DSepDeterminer, CACalculator, varargin)
%GLGLEARNER Create max-min local-to-global learner.
%   L = ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.MMLGLEARNER(D, C)
%   creates a max-min local-to-global learner with d-separation determiner
%   D and conditional-association calculator C.
%
%   [...] = ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.MMLGLEARNER(D, C,
%   'Param1', VAL1, 'Param2', VAL2, ...) specifies additional parameter
%   name/value pairs as in
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.GLGLEARNER/GLGLEARNER.
%
%   See also ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.DSEPDETERMINER,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.MM.CACALCULATOR,
%   ORG.MENSXMACHINA.PGM.BN.LEARNING.CB.LG.GL.GLGLEARNER/GLGLEARNER.
    
    % call GLGLEARNER constructor
    Obj = Obj@glglearner(DSepDeterminer, varargin{:});
    
    % parse specific input
    assert(isa(CACalculator, 'cacalculator') && isscalar(CACalculator) && CACalculator.nVars == DSepDeterminer.nVars);
    
    % set specific properties
    Obj.CACalculator = CACalculator;

end

end

methods(Access = protected)

% abstract method implementations

function xTpc = initializetpc(~, ~)
    
    % initialize empty TPC(X)
    xTpc = zeros(0, 1);
    
end

function xOpenPriorityQueue = inititializeopenpriorityqueue(~, xOpenInd)

    % initialize minimal associations with maximal values
    xOpenPriorityQueue = struct(...
        'minCA', Inf(size(xOpenInd)), ...
        'minCACA2', Inf(size(xOpenInd)));

end

function [xOpenInd xOpenPriorityQueue] = ...
        updateopen(Obj, xInd, xTpcInd, xOpenInd, xOpenPriorityQueue)

    import sepsetidentifieddata;
    import mmlglearner;

    % calculates associations of the remaining variables in OPEN with Y 
    % given the new subsets of TPC(Y), while returning any variables
    % that achieve no association
    
    % initialize members of OPEN(x) that have zero CA
    xOpenMemberHasZeroCA = false(size(xOpenInd));

    for yIndInXOpen = 1:length(xOpenInd) % for each variable Y in OPEN(X)

        % get index in V
        yInd = xOpenInd(yIndInXOpen);

        % upper-bound sepset cardinality
        msc = maxsepsetcard(Obj.DSepDeterminer, xInd, yInd, xTpcInd);
        
        if isnan(msc)
            msc = 0; % allow CA(Y;X) to be calculated
        end

        % limit maximal sepset cardinality
        msc = min(msc, Obj.maxSepsetCard);

        if isempty(xTpcInd)
            
            % calculate CA(Y;X)
            [thisCA thisCA2] = ca(Obj.CACalculator, xInd, yInd, zeros(1, 0));
                
%             fprintf('\nCalculated CA(%d,%d) = %g\n', xInd, yInd, thisCA);

            if thisCA == 0

                xOpenMemberHasZeroCA(yIndInXOpen) = true;
                    
                notify(Obj, 'sepsetIdentified', sepsetidentifieddata(xInd, yInd, zeros(1, 0)));
                
            elseif mmlglearner.compareca(thisCA, thisCA2, xOpenPriorityQueue.minCA(yIndInXOpen), xOpenPriorityQueue.minCACA2(yIndInXOpen)) < 0

                % update minimum CA and corresponding secondary CA
                xOpenPriorityQueue.minCA(yIndInXOpen) = thisCA;
                xOpenPriorityQueue.minCACA2(yIndInXOpen) = thisCA2;

            end

            continue;

        end

%         fprintf('\nmax |Z| for a subset Z of TPC(%d) = [%s] in testing X(%d;%d|Z): %d of %d\n', xInd, num2str(xTpcInd), xInd, yInd, msc, length(xTpcInd));

        if msc == 0 

            % only X(Y;X|[]) can be tested accuratelly
            % and it's already tested
            continue;

        end % else msc > 0
        
        % get F index
        fInd = xTpcInd(end);

        % calculate CA(Y;X|F)
        [thisCA thisCA2] = ca(Obj.CACalculator, xInd, yInd, fInd);
                
%         fprintf('\nCalculated CA(%d,%d|[%d]) = %g\n', xInd, yInd, fInd, thisCA);

        if thisCA == 0
            
            xOpenMemberHasZeroCA(yIndInXOpen) = true;
                    
            notify(Obj, 'sepsetIdentified', sepsetidentifieddata(xInd, yInd, fInd));
            
            continue; % no need to try more subsets
            
        elseif mmlglearner.compareca(thisCA, thisCA2, xOpenPriorityQueue.minCA(yIndInXOpen), xOpenPriorityQueue.minCACA2(yIndInXOpen)) < 0

            % update minimum association and corresponding secondary association
            xOpenPriorityQueue.minCA(yIndInXOpen) = thisCA;
            xOpenPriorityQueue.minCACA2(yIndInXOpen) = thisCA2;

        end

        if msc < 2
            % |Z| == 0 ([]) and |Z| == 1 (F) already tested,
            % no other |Z| < 2 to test
            continue;                
        end

        for zCard = 2:msc % for each cardinality up to the maximum

            % find subsets of this size - 1 among the previous TPC(Y)
            C = nchoosek( xTpcInd(1:end-1), zCard - 1); 

            % note: when |TPC(Y)| == 2, msc == 2 (due to the condition above),
            % length(xTpcInd(1:end-1)) == 1 and zCard - 1 == 1,
            % nchoosek( xTpcInd(1:end-1), zCard - 1) == nchoosek( xTpcInd(1:end-1), 1)
            % == xTpcInd(1:end-1)

            % debug
    %                     if length(xTpcInd) == 2
    %                         assert(msc == 2);
    %                         assert(length(xTpcInd(1:end-1)) == 1);
    %                         assert(zCard - 1 == 1);
    %                         assert(nchoosek( xTpcInd(1:end-1), zCard - 1) == xTpcInd(1:end-1));
    %                     end        

            for k = 1:size(C, 1) % for each such subset

                % add F and let it be Z 
                zInd = [ C(k, :) xTpcInd(end)];

                % xInd.e. for each subset Z of TPC(Y) including F
                % and with cardinality up to the maximum

                % calculate CA(Y;X|Z)                
                [thisCA thisCA2] = ca(Obj.CACalculator, xInd, yInd, zInd);
                
%                 fprintf('\nCalculated CA(%d,%d|[%s]) = %g\n', xInd, yInd, num2str(zInd), thisCA);

                if thisCA == 0
                    
                    xOpenMemberHasZeroCA(yIndInXOpen) = true;
                    
                    notify(Obj, 'sepsetIdentified', sepsetidentifieddata(xInd, yInd, zInd));

                    break; % no need to try more subsets

                elseif mmlglearner.compareca(thisCA, thisCA2, xOpenPriorityQueue.minCA(yIndInXOpen), xOpenPriorityQueue.minCACA2(yIndInXOpen)) < 0

                    % update minimum association and corresponding secondary association
                    xOpenPriorityQueue.minCA(yIndInXOpen) = thisCA;
                    xOpenPriorityQueue.minCACA2(yIndInXOpen) = thisCA2;

                end

            end

            if xOpenMemberHasZeroCA(yIndInXOpen)
                break;
            end

        end

    end
    
%     fprintf('\nRemoved [%s] from OPEN(%d)\n', num2str(xOpenInd(xOpenMemberHasZeroCA)), xInd);

    % throw variables with no association with X given TPC(X) from OPEN    
    xOpenInd(xOpenMemberHasZeroCA) = [];
    xOpenPriorityQueue.minCA(xOpenMemberHasZeroCA) = [];
    xOpenPriorityQueue.minCACA2(xOpenMemberHasZeroCA) = [];

end

function [xTpcInd xOpenInd xOpenPriorityQueue] = ...
    removefromopen(~, ~, xTpcInd, xOpenInd, xOpenPriorityQueue)

    import mmlglearner;

    % Max-min heuristic

    if isempty(xOpenInd)

        % nothing to remove
        return;

    end

    % get the index (among the remaining variables) and minimum association
    % of the variable with the maximum minimum association
    % among the remaining variables

    maxMinCA = 0;
    maxMinCACA2 = 0;

    for yIndInXOpen = 1:numel(xOpenInd) % for each variable Y in OPEN(X)

        if mmlglearner.compareca(xOpenPriorityQueue.minCA(yIndInXOpen), xOpenPriorityQueue.minCACA2(yIndInXOpen), maxMinCA, maxMinCACA2) > 0

            fIndInXOpen = yIndInXOpen;
            
            maxMinCA = xOpenPriorityQueue.minCA(yIndInXOpen);
            maxMinCACA2 = xOpenPriorityQueue.minCACA2(yIndInXOpen);

        end

    end

    % insert it in TPC(X)
    xTpcInd = [xTpcInd xOpenInd(fIndInXOpen)];
                    
%     fprintf('\nMoved [%s] from OPEN(%d) to TPC(%d)\n', num2str(xOpenInd(fIndInXOpen)), xInd, xInd);

    % remove it from OPEN(X)
    xOpenInd(fIndInXOpen) = [];
    xOpenPriorityQueue.minCA(fIndInXOpen) = [];
    xOpenPriorityQueue.minCACA2(fIndInXOpen) = [];

end

function xTpcInd = removefromtpc(Obj, xInd, xTpcInd, xOpenInd)

    %PHASE2 Phase II (backward): remove variables from TPC(X)

    import sepsetidentifieddata;
    
    if ~isempty(xOpenInd)
        
        % OPEN(X) not empty yet...
        return;

    end

    if length(xTpcInd) < 2

        % if 0, nothing to remove
        % if 1, there are no other variables in TPC(X) that could be affected
        % by the insertion of the last and only variable F
        return;

    end
    
    xTpcMemberIsRemoved = false(size(xTpcInd));

    for yIndInXTpc = 1 : length(xTpcInd) - 1 % for each variable Y in TPC(X) except the last one (F)

        % let it be Y
        yInd = xTpcInd(yIndInXTpc);

        % we exclude F, the last variable added to TPC(X),
        % because it is already known that F is associated with X
        % given any subset of the previous TPC(X)

        % find indices of the rest variables in TPC(X) not yet removed;
        % in the extreme case this in only the last variable (F)

        retainedXTpcButYIndInXTpc = setdiff([ 1 : (yIndInXTpc - 1) (yIndInXTpc + 1) : length(xTpcInd)], find(xTpcMemberIsRemoved));
        
        % find the maximal conditioning set cardinality such that the test is reliable
        msc = maxsepsetcard(Obj.DSepDeterminer, xInd, yInd, xTpcInd(retainedXTpcButYIndInXTpc));
        
%         assert(~isnan(msc));
   
        % get minimum between reliability criterion maxk and user MAXK
        msc = min(msc, Obj.maxSepsetCard);
        
        if msc == 0
            continue;
        end

        for zCard = 1:msc % for each cardinality up to the maximum

            % note: for Y in TPC(X), X(Y;X|[]) has been tested; 
            % ignoring size 0

            % find subsets of -indices- of variables in TPC(X) 
            % of this size that do not contain the index of Y      
            C = nchoosek(retainedXTpcButYIndInXTpc, zCard);

            % note: when length(indices) == 1, msc == 1 (due to the condition
            % above), and zCard == 1, nchoosek(indices, zCard)
            % == nchoosek(indices, 1) == indices

%                 % debug
%                 if length(indices) == 1
%                     assert(msc == 1);
%                     assert(zCard == 1);
%                     assert(nchoosek(indices, zCard) == indices);
%                 end

            for k = 1:size(C, 1) % for each such index subset

                if max(C(k, :)) < yIndInXTpc

                    % subset is older in TPC(X) than Y, test already performed
                    continue;

                end

                % let it be Z
                zInd = xTpcInd(C(k, :));
                
%                 fprintf('\nChecking CI(%d,%d|%s)...', xInd, yInd, num2str(zInd));
                
                if isdsep(Obj.DSepDeterminer, xInd, yInd, zInd)

                    xTpcMemberIsRemoved(yIndInXTpc) = true;
                    
                    notify(Obj, 'sepsetIdentified', sepsetidentifieddata(xInd, yInd, zInd));
                
%                     fprintf(' Yes.\n');

                    break; % no need to try more subsets
                    
                end
                
%                 fprintf(' No.\n');

            end

            if xTpcMemberIsRemoved(yIndInXTpc)
                break;
            end

        end

    end
    
%     fprintf('\nRemoved [%s] from TPC(%d)\n', num2str(xTpcInd(xTpcMemberIsRemoved)), xInd);
                
    % remove variables from TPC(X)
    xTpcInd(xTpcMemberIsRemoved) = [];

end

function tf = mustterminate(~, ~, ~, xOpenInd)

    tf = isempty(xOpenInd);

end

end

methods(Static, Access = private)

function cmp = compareca(ca1, ca21, ca2, ca22)

    % compare conditional associations

    if ca1 < ca2
        cmp = -1; 
    elseif ca1 > ca2
        cmp = 1;
    else

        if ca21 < ca22
            cmp = -1;
        elseif ca21 > ca22
            cmp = 1;
        else
            cmp = 0;
        end

    end

end

end

end