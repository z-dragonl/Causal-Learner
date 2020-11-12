classdef(Sealed) heuristicapplier < citrcapplier
%HEURISTICAPPLIER Heuristic-power-rule applier.
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.HEURISTIC.HEURISTICAPPLIER is the
%   class of heuristic-power-rule appliers. The heuristic power rule is a
%   reliability criterion for Chi-square tests of conditional independence.
%   According to this criterion, a test is reliable if there are at least
%   h-ps observations, on average, in each cell of the contigency table of
%   the test.

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

% References:
% [1] Aliferis, C. F., Statnikov, A., Tsamardinos, I., Mani, S., and
%     Koutsoukos, X. D. 2010. Local Causal and Markov Blanket Induction for
%     Causal Discovery and Feature Selection for Classification Part I:
%     Algorithms and Empirical Evaluation. J. Mach. Learn. Res. 11 (Mar.
%     2010), 171-234.

properties(SetAccess = immutable)
    
nVars % number of variables -- a numeric positive integer
    
end

properties(GetAccess = private, SetAccess = immutable)
    
nObs % number of observations -- a numeric positive integer

varNValues % variable numbers of values -- an 1-by-N numeric vector of positive integers 

hps % h-ps -- a numeric nonnegative integer

end

methods

% constructor

function Obj = heuristicapplier(sample, varNValues, hps)
%HEURISTICAPPLIER Create heuristic-power-rule applier.
%   OBJ =
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.HEURISTIC.HEURISTICAPPLIER(SAMPLE,
%   VARNVALUES) creates a heuristic-power-rule applier with sample SAMPLE,
%   variable numbers of values VARNVALUES and h-ps 5. SAMPLE is an M-by-N
%   numeric matrix of positive integers, where M is the number of
%   observations and N is the number of variables in the sample. Each
%   element of SAMPLE is the linear index of the value of the corresponding
%   variable for the corresponding observation. VARNVALUES is an 1-by-N
%   numeric vector of positive integers. Each element of VARNVALUES is the
%   number of values of the corresponding variable.
%
%   OBJ =
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.HEURISTIC.HEURISTICAPPLIER(SAMPLE,
%   VARNVALUES, HPS) sets h-ps to HPS. HPS is a numeric nonnegative
%   integer.

    import org.mensxmachina.stats.array.*;

    % parse input
    validateattributes(sample, {'numeric'}, {'2d', 'positive', 'integer'});
    validateattributes(varNValues, {'numeric'}, {'size', [1 size(sample, 2)], 'positive', 'integer'});    
    
    if nargin < 3
        hps = 5;
    else
        validateattributes(hps, {'numeric'}, {'scalar', 'nonnegative', 'integer'});
    end
    
    % set properties
    Obj.nVars = size(sample, 2);
    Obj.nObs = size(sample, 1);
    Obj.varNValues = varNValues;
    Obj.hps = hps;

end

% abstract method implementations

function tf = isreliablecit(Obj, i, j, k)
    
    % (no validation)
    
    % get test variable numbers of values
    testVarNValues = Obj.varNValues([i j k]);
    
    tf = Obj.nObs / prod(testVarNValues) >= Obj.hps;
    
end

function wmk = worstmaxcondsetcard(Obj, i, j, k)
    
    % (no validation)
    
    wmk = maxcondsetcard(Obj, i, j, k, 'descend');
    
end

function bmx = bestmaxcondsetcard(Obj, i, j, k)
    
    % (no validation)
    
    bmx = maxcondsetcard(Obj, i, j, k, 'ascend');
    
end

end

methods(Access = private)
    
function maxcondsetcard = maxcondsetcard(Obj, i, j, k, sNValuesSortMode)
    
    % get number of values for each test variable
    xNValues = Obj.varNValues(i);
    yNValues = Obj.varNValues(j);
    sNValues = Obj.varNValues(k);    
    
    testVarNValues = [xNValues yNValues sort(sNValues, sNValuesSortMode)];

%     if length(testVarNValues) == 2
%         maxcondsetcard = 0;
%         return;
%     end
% 
%     % conservative estimate of MAXK
%     % derived from nobs/(testVarNValues(j1)*testVarNValues(j2)*min(testVarNValues(j3))^maxcondsetcard) >= hps
%     maxcondsetcard = round(log(Obj.nObs/(Obj.hps*testVarNValues(1)*testVarNValues(2)))/log(min(testVarNValues(3:end))));
% 
%     % limit maxcondsetcard to [0, length(j3)]
%     maxcondsetcard = max(min(maxcondsetcard, length(testVarNValues) - 2), 0);
% 
%     return;

    maxmaxk = length(testVarNValues) - 2;

    % N/PROD(NUMLEVELS(1:(2 + K))) >= HPS solved for K

    maxcondsetcard = 0;

    while maxcondsetcard <= maxmaxk && Obj.nObs/prod(testVarNValues(1:(2 + maxcondsetcard))) >= Obj.hps
        maxcondsetcard = maxcondsetcard + 1;
    end

    maxcondsetcard = maxcondsetcard - 1;

    if maxcondsetcard < 0
        maxcondsetcard = NaN;
    end

end

end

end