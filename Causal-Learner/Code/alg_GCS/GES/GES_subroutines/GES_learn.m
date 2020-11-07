
% Greedy equivalence search
function [Record] = GES_learn(X,score_type,maxP,parameters)
% INPUT:
% X: Data with T*D dimensions
% score_type: the score function you want to use
% maxP: allowed maximum number of parents when searching the graph
% parameters: when using CV likelihood, 
%               parameters.kfold: k-fold cross validation
%               parameters.lambda: regularization parameter

% OUTPUT:
% Record.G: learned causal graph
% Record.update1: each update (Insert operator) in the forward step
% Record.update2: each update (Delete operator) in the backward step
% Record.G_step1: learned graph at each step in the forward step
% Record.G_step2: learned graph at each step in the backward step
% Record.score: the score of the learned graph


if(score_type == 1) % k-fold negative cross validated likelihood based on regression in RKHS
    score_func = 'local_score_CV_general'; 
    if(nargin<4)
        parameters.kfold = 10; % 10 fold cross validation
        parameters.lambda = 0.01; % regularization parameter
    end
end
if(score_type == 2) % negative marginal likelihood based on regression in RKHS
    score_func = 'local_score_marginal_general';
    parameters = [];
end

if(nargin<3)
    maxP = size(X,2)/2; % maximum number of parents
end

N = size(X,2); % number of variables
G = zeros(N,N); % initialize the graph structure
score = Score_G(X,G,score_func,parameters); % initialize the score
G = PDAG2DAG(G);
G = DAG2CPDAG(G); % transform the input DAG to a CPDAG

%%  forward greedy search
record_local_score=cell(1,N); % record the local score calculated each time. Thus when we transition to the second phase, many of the operators can be scored without an explicit call the the scoring function
%record_local_score{trial}{j} record the local scores when Xj as a parent
score_new=score;
count1=0;
update1=[];
G_step1=[];
while(1)
    count1=count1+1;
    score=score_new;
    score_record1(count1)=score;
    graph_record1{count1}=G;
    min_chscore = 1e7;
    min_desc={};
    for i=1:N
        for j=1:N
            if (G(i,j)==0 & G(j,i)==0 & i~=j & length(find(G(:,j)==1))<=maxP)  % find a pair (Xi, Xj) that is not adjacent in the current graph , and restrict the number of parents
                Tj = find(G(j,:)==-1); % neighbors of Xj
                Ti = union(find(G(i,:)~=0),find(G(:,i)~=0)'); % adjacent to Xi
                NTi = setdiff(1:N,Ti);
                T0 = intersect(Tj,NTi); % find the neighbours of Xj that are not adjacent to Xi
                % for any subset of T0
                [sub] = Combinatorial (T0); % find all the subsets for T0
                S = zeros(1,length(sub));
                % S indicate whether we need to check sub{k}.
                % 0: check both conditions.
                % 1: only check the first condition
                % 2: check nothing and is not valid.
                for k=1:length(sub)
                    if(S(k)<2) % S indicate whether we need to check subset(k)
                        V1 = Insert_validity_test1(G,i,j,sub{k}); % Insert operator validation test:condition 1
                        if(V1)
                            if(~S(k))
                                V2 = Insert_validity_test2(G,i,j,sub{k});% Insert operator validation test:condition 2
                            else
                                V2=1;
                            end
                            if(V2)
                                Idx = find_subset_include(sub{k},sub);% find those subsets that include sub(k)
                                S(find(Idx==1))=1;
                                
                                [chscore,desc,record_local_score] = Insert_changed_score(X,G,i,j,sub{k},record_local_score,score_func,parameters);%calculate the changed score after Insert operator
                                % desc{count} saves the corresponding (i,j,sub{k})
                                if(chscore<min_chscore)
                                    min_chscore = chscore;
                                    min_desc = desc;
                                end
                            end
                        else
                            Idx = find_subset_include(sub{k},sub);% find those subsets that include sub(k)
                            S(find(Idx==1))=2;
                        end
                    end
                end
            end
        end
    end
    
    if(~isempty(min_desc))
        score_new = score+min_chscore;
        if(score-score_new <= 0)
            break;
        end
        G = Insert(G,min_desc{1},min_desc{2},min_desc{3});
        update1{count1}=[min_desc{1},min_desc{2},min_desc{3}];
        G = PDAG2DAG(G);
        G = DAG2CPDAG(G);
        G_step1{count1}=G;
    else
        score_new = score;
        break;
    end
end


%% backward greedy search
count2=0;
score_new=score;
update2=[];
G_step2=[];
while(1)
    count2=count2+1;
    score=score_new;
    score_record2(count2)=score;
    graph_record2{count2}=G;
    min_chscore = 1e7;
    min_desc={};
    for i=1:N
        for j=1:N
            if((G(i,j)==-1 | G(i,j)==1)) % if Xi - Xj or Xi -> Xj
                Hj = find(G(j,:)==-1); % neighbors of Xj
                Hi = union(find(G(i,:)~=0),find(G(:,i)~=0)'); % adjacent to Xi
                H0 = intersect(Hj,Hi); % find the neighbours of Xj that are adjacent to Xi
                % for any subset of H0
                [sub] = Combinatorial (H0); % find all the subsets for H0
                S = ones(1,length(sub)); % S indicate whether we need to check sub{k}.
                % 1: check the condition,
                % 2: check nothing and is valid;
                for k=1:length(sub)
                    if(S(k)==1)
                        V = Delete_validity_test(G,i,j,sub{k}); % Delete operator validation test
                        if(V)
                            % find those subsets that include sub(k)
                            Idx = find_subset_include(sub{k},sub);
                            S(find(Idx==1))=2; %and set their S to 2
                        end
                    else
                        V=1;
                    end
                    if(V)
                        [chscore,desc,record_local_score] = Delete_changed_score(X,G,i,j,sub{k},record_local_score,score_func,parameters); %calculate the changed score after Insert operator
                        % desc{count} saves the corresponding (i,j,sub{k})
                        if(chscore<min_chscore)
                            min_chscore = chscore;
                            min_desc = desc;
                        end
                    end
                end
            end
        end
    end
    
    if(~isempty(min_desc))
        score_new = score+min_chscore;
        if(score-score_new <=0 )
            break;
        end
        G = Delete(G,min_desc{1},min_desc{2},min_desc{3});
        update2{count2} = [min_desc{1},min_desc{2},min_desc{3}];
        G = PDAG2DAG(G);
        G = DAG2CPDAG(G);
        G_step2{count2}=G;
    else
        score_new = score;
        break;
    end
end
Record.update1=update1;
Record.update2=update2;
Record.G_step1 = G_step1;
Record.G_step2 = G_step2;
Record.G=G;
Record.score = score;
end
