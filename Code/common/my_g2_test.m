function [pvalue, dep] = my_g2_test(X, Y, S, Data, ns, alpha)


[nSamples,nVars]=size(Data);

if isreliablecit(ns, nSamples, X, Y, S) % test is reliable, attempt it
    
    % calculate test p-value
    [p,stat] = citpvalue(Data,ns, X, Y, S);
    
    % compare p-value with alpha
    
    if p > alpha % independent
        
        ca = 0;
        ca2 = 0;
        
    else % dependent
        
        ca = 1/p;
        ca2 = abs(stat);
        
    end
    
else
    
    p = NaN;
    stat = NaN;
    
    if isempty(S)
        ca = 0;
        ca2 = 0;
    else
        ca = Inf;
        ca2 = Inf;
    end
    

end


pvalue=p;
dep=ca2;

end



function [p,stat] = citpvalue(Data, ns, i, j, k)

% (no validation)

ind = [i j k];

% select test variables
testSample = Data(:, ind);
testNLevels = ns(:, ind);

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
stat = chi2stat(Obs_vector, Exp_vector);

% call gammainc to avoid roundoff in computing upper tail probability
p = gammainc(stat/2, df/2, 'upper');  % p = 1 - chi2cdf(stat,df);

end


function stat = chi2stat(obs, exp)

% compute G statistic

terms = obs.*log(obs./exp);

% set NaN terms to zero
terms(isnan(terms)) = 0;

stat = 2*sum(terms);

end



function tf = isreliablecit(ns, nSamples,i, j, k)

hps=5; %default

% (no validation)

% get test variable numbers of values
testVarNValues = ns([i j k]);

tf = nSamples / prod(testVarNValues) >= hps;

end
