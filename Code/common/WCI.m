function [pvalue, dep,stat] = WCI(X, Y, S, Data, ns, alpha)

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

% i=X;
% j=Y;
% k=S;

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
stat = chi2stat_w(Obs_vector, Exp_vector,df, Data);

% call gammainc to avoid roundoff in computing upper tail probability
p = gammainc(stat/2, df/2, 'upper');  % p = 1 - chi2cdf(stat,df);

end

%  obs=Obs_vector;
%  exp=Exp_vector;


function stat = chi2stat_w(obs, exp,df, Data)


[~,nVars]=size(Data);

obs_sum=sum(obs);
obs_ce=obs./obs_sum;
ce=-log(obs_ce);
ce(obs==0)=0;
not_z_ce=ce(ce~=0);
exp_ce=mean(not_z_ce);
cs_ce=(ce-exp_ce)./exp_ce;
% cs_ce(cs_ce>0)=-cs_ce(cs_ce>0);
cs_ce(obs==0)=0;
a=0.01;


if  nVars == 109
    
terms = obs.*log(obs./exp); % 
% set NaN terms to zero
terms(isnan(terms)) = 0;
terms=terms.*2;
temp_stat = sum(terms);
temp_p = gammainc(temp_stat/2, df/2, 'upper');

else
    
terms = obs.*log(obs./exp); % 
% set NaN terms to zero
terms(isnan(terms)) = 0;
terms=terms.*2;
temp_stat = sum(terms);
temp_p=0;

end



x_terms = ((obs-exp).^2)./exp;
y_terms = ((abs(obs-exp)-0.5).^2)./exp;
noise=abs(x_terms-y_terms);
noise(obs<exp)=-noise(obs<exp);
noise(obs==0)=0;
exp_zore=exp(obs==0);
exp_zz= (exp_zore~=0);


if temp_p > a % independent 

     abs_obs_exp=abs(obs-exp);
     abs_obs_exp_ln=~(abs_obs_exp>=0.5);
     terms(abs_obs_exp_ln)=terms(abs_obs_exp_ln)+noise(abs_obs_exp_ln);
     abs_obs_exp_gn=(abs_obs_exp>=0.5);
     terms(abs_obs_exp_gn)=terms(abs_obs_exp_gn)-noise((abs_obs_exp_gn));
     stat = sum(terms)

else % dependent 

      if ~isempty( exp_zz)
          terms = terms+terms.*cs_ce;
      end
       stat = sum(terms)*0.94;

end  % end p>alpha 


end



function tf = isreliablecit(ns, nSamples,i, j, k)

hps=5; %default  
% (no validation)
% get test variable numbers of values
testVarNValues = ns([i j k]);
tf = nSamples / prod(testVarNValues) >= hps;

end
