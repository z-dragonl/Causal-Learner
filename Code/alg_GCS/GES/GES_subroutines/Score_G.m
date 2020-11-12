

function score = Score_G(Data, G, score_func, parameters)  % calculate the score for the current G
% here G is a DAG
N=size(G,1);
score = 0;
for i=1:N
    PA = find(G(:,i)==1)';
    delta_score = feval(score_func,Data,i,PA,parameters);
    score = score + delta_score; % also need to consider the case when PA is empty!!
end
