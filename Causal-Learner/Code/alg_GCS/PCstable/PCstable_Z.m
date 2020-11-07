function [DAG,time,test] = PCstable_Z(Data,alpha,samples,p,maxK)
%
% PCstable_Z learns a causal graph on continuous data
%
% INPUT :
%       Data is the data matrix
%       alpha is the significance level
%       samples is the number of data samples
%       p is the number of nodes
%       maxK is the maximum size of conditioning set
%
% OUTPUT:
%       DAG is the causal graph
%       time is the runtime of the algorithm
%       test is the number of conditional independence tests
%
%


if (nargin == 2)
   [samples,p]=size(Data);
   maxK=3;
end

start=tic;
sep = cell(p,p);
ord = 0;
done = 0;
G = ones(p,p);

G=setdiag(G,0);
test=0;


while ~done && ord<=maxK
    done = 1;
    [X,Y] = find(G);
    
    % LINES 5-7 of pseudo code
    ADJ=cell(1,length(X));
    for i=1:length(X)
        x = X(i); y = Y(i);
        nbrs = mysetdiff(myneighbors(G, y), x);    % bug fix by Raanan Yehezkel <raanany@ee.bgu.ac.il> 6/27/04
        ADJ{i}=nbrs;
    end
    % LINES 5-7 of pseudo code   END

    for i=1:length(X)
        x = X(i); y = Y(i);

        if length(ADJ{i}) >= ord && G(x,y) ~= 0 
            done = 0;

            SS = subsets1(ADJ{i}, ord);
            for si=1:length(SS)
                S = SS{si};
                test=test+1;
                
                [CI]=my_fisherz_test(x,y,S,Data,samples,alpha);
                if isnan(CI)
                    CI=0;
                end
                
                if(CI==1)

                    G(x,y) = 0;
                    G(y,x) = 0;

                    sep{x,y} = myunion(sep{x,y}, S);
                    sep{y,x} = myunion(sep{y,x}, S);
                    break; % no need to check any more subsets
                end
            end
        end
    end
    ord = ord + 1;
end


% disp('finishing skeletons,creating V structures....');

% Create the minimal pattern,
% i.e., the only directed edges are V structures.
DAG = G;                 
[X, Y] = find(G);
% We want to generate all unique triples x,y,z
% This code generates x,y,z and z,y,x.
for i=1:length(X)
    x = X(i);
    y = Y(i);
    Z = find(G(y,:));
    Z = mysetdiff(Z, x);
    for z=Z(:)'
        if G(x,z)==0 && ~ismember(y, sep{x,z}) && ~ismember(y, sep{z,x})
            %fprintf('%d -> %d <- %d\n', x, y, z);
            DAG(x,y) = -1; DAG(y,x) = 0;
            DAG(z,y) = -1; DAG(y,z) = 0;
        end
    end
end

% disp('finishing V structures, directed edge directions....');

% Convert the minimal pattern to a complete one,
% i.e., every directed edge in P is compelled
% (must be directed in all Markov equivalent models),
% and every undirected edge in P is reversible.
% We use the rules of Pearl (2000) p51 (derived in Meek (1995))

old_pdag = zeros(p);
iter = 0;
while ~isequal(DAG, old_pdag)
    iter = iter + 1;
    old_pdag = DAG;
    % rule 1            % a -> b --C  => b->c
    [A,B] = find(DAG==-1); % a -> b
    for i=1:length(A)
        a = A(i); b = B(i);
        C = find(DAG(b,:)==1 & G(a,:)==0); % all nodes adj to b but not a
        if ~isempty(C)
            DAG(b,C) = -1; DAG(C,b) = 0;
            %fprintf('rule 1: a=%d->b=%d and b=%d-c=%d implies %d->%d\n', a, b, b, C, b, C);
        end
    end
    % rule 2            % a->c->b, a--b => a->b
    [A,B] = find(DAG==1); % unoriented a-b edge
    for i=1:length(A)
        a = A(i); b = B(i);
        if any( (DAG(a,:)==-1) & (DAG(:,b)==-1)' );
            DAG(a,b) = -1; DAG(b,a) = 0;
            %fprintf('rule 2: %d -> %d\n', a, b);
        end
    end
    % rule 3            % a--c->b, a--d->b, pdag(c,d)=pdag(d,c)=0, a--b  => a->b
    [A,B] = find(DAG==1); % a-b
    for i=1:length(A)
        a = A(i); b = B(i);
        C = find( (DAG(a,:)==1) & (DAG(:,b)==-1)' );
        % C contains nodes c s.t. a-c->ba
        G2 = setdiag(G(C, C), 1);
        if any(G2(:)==0) % there are 2 different non adjacent elements of C
            DAG(a,b) = -1; DAG(b,a) = 0;
            %fprintf('rule 3: %d -> %d\n', a, b);
        end
    end
end


DAG(DAG==-1)=1;

DAG=cpdag_to_dag(DAG);

% draw_graph(DAG);

time=toc(start);


