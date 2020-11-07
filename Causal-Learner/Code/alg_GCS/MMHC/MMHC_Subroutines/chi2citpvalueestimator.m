classdef chi2citpvalueestimator < citpvalueestimator
%CHI2CITPVALUEESTIMATOR Chi-square-test-of-conditional-independence-p-value estimator.
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.CHI2CITPVALUEESTIMATOR is the
%   abstract class of Chi-square-test-of-conditional-independence-p-value
%   estimators. The statistics used in Chi-square tests of conditional
%   independence follow the Chi-square distribution.
%
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.CHI2CITPVALUEESTIMATOR/CITPVALUE
%   returns P = 1 and STAT = 0 when the degrees of freedom of the test are
%   zero.

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
%   [1] http://en.wikipedia.org/wiki/Chi-squared_test

properties(SetAccess = immutable)
    
nVars % number of variables -- a numeric nonnegative integer
    
end

properties(GetAccess = private, SetAccess = immutable)
    
sample % sample -- an M-by-N numeric matrix of positive integers

varNValues % variable numbers of values -- an 1-by-N numeric vector of positive integers

end

methods

% constructor

function Obj = chi2citpvalueestimator(sample, varNValues)
%CHI2CITPVALUEESTIMATOR Create Chi-square-test-of-conditional-independence-p-value estimator.
%   OBJ =
%   ORG.MENSXMACHINA.STATS.TESTS.CI.CHI2.CHI2CITPVALUEESTIMATOR(SAMPLE,
%   VARNVALUES) creates a Chi-square-test-of-conditional-independence
%   p-value estimator with sample SAMPLE and variable numbers of values
%   VARNVALUES. SAMPLE is an M-by-N numeric matrix of positive integers,
%   where M is the number of observations and N is the number of variables
%   in the sample. Each element of SAMPLE is the linear index of the value
%   of the corresponding variable for the corresponding observation.
%   VARNVALUES is an 1-by-N numeric vector of positive integers. Each
%   element of VARNVALUES is the number of values of the corresponding
%   variable.


    % parse input
    validateattributes(sample, {'numeric'}, {'2d', 'positive', 'integer'});
    validateattributes(varNValues, {'numeric'}, {'size', [1 size(sample, 2)], 'positive', 'integer'});
    
    % set properties
    Obj.nVars = size(sample, 2);
    Obj.sample = sample;
    Obj.varNValues = varNValues;

end

end

methods(Sealed)

% abstract method implementations

function [p stat] = citpvalue(Obj, i, j, k)
    
    % (no validation)

    ind = [i j k];
    
    % select test variables
    testSample = Obj.sample(:, ind);
    testNLevels = Obj.varNValues(:, ind);

    nObs = size(testSample, 1);
    nTestVars = size(testSample, 2);
    
    % compute observed level combination counts
    Obs = accumarray(testSample, ones(1, nObs), testNLevels);

    % Obs_xs(i,j,kk{:}): observed count of (i,kk{:}), same for all j
    ObsSum2 = sum(Obs, 2); % sum for all levels of j
    Obs_xs_size = ones(1, nTestVars);
    Obs_xs_size(2) = testNLevels(2);
    Obs_xs = repmat(ObsSum2, Obs_xs_size); % replicate for all levels of j

    % Obs_ys(i,j,kk{:}): observed count of (j,kk{:}), same for all i
    ObsSum1 = sum(Obs, 1); % sum for all levels of i
    Obs_ys_size = ones(1, nTestVars);
    Obs_ys_size(1) = testNLevels(1);
    Obs_ys = repmat(ObsSum1, Obs_ys_size); % replicate for all levels of i

    % Obs_s(i,j,kk{:}): observed count of (kk{:}), same for all i and j
    Obs_s = sum(ObsSum1, 2); % sum for all levels of i and j
    Obs_s_size = ones(1, nTestVars);
    Obs_s_size([1 2]) = testNLevels([1 2]);
    Obs_s = repmat(Obs_s, Obs_s_size); % replicate for all levels of i and j

    % compute expected level combination counts
    Exp = Obs_xs.*Obs_ys./Obs_s;

    j3NLevelsProd = prod(testNLevels(3:end));
    
    ObsSum1 = reshape(ObsSum1, testNLevels(2), j3NLevelsProd);
    ObsSum2 = reshape(ObsSum2, testNLevels(1), j3NLevelsProd);

    % note: for seems faster
    df = 0;

    for iComb = 1:j3NLevelsProd

        df = df + max(testNLevels(1) - 1 - sum(~ObsSum2(:,iComb)), 0) * max(testNLevels(2) - 1 - sum(~ObsSum1(:,iComb)), 0);
    
    end

    if df == 0
        
        % independence
        p = 1;
        stat = 0;

        return;

    end

    Obs_vector = Obs(:);
    Exp_vector = Exp(:);

    % compute test statistic
    stat = Obj.chi2stat(Obs_vector, Exp_vector);

    % call gammainc to avoid roundoff in computing upper tail probability
    p = gammainc(stat/2, df/2, 'upper');  % p = 1 - chi2cdf(stat,df);

end

end

methods(Abstract, Access = protected)
 
%CHI2STAT Compute Chi-square statistic.
stat = chi2stat(obs, exp);

end

end