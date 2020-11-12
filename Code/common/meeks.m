
function   [DAG,pdag,G]=meeks(DAG,pdag,G,p)


old_pdag = zeros(p);
iter = 0;

while ~isequal(pdag, old_pdag)
    iter = iter + 1;
    old_pdag = pdag;
    % rule 1
    [A,B] = find(pdag==-1); % a -> b
    for i=1:length(A)
        a = A(i); b = B(i);
        C = find(pdag(b,:)==1 & DAG(a,:)==0); % all nodes adj to b but not a
        if ~isempty(C)       
            pdag(b,C) = -1; pdag(C,b) = 0;
            G(b,C) = 1; G(C,b) = 0;
            %fprintf('rule 1: a=%d->b=%d and b=%d-c=%d implies %d->%d\n', a, b, b, C, b, C);
        end
    end
    

    % rule 2
    [A,B] = find(pdag==1); % unoriented a-b edge
    for i=1:length(A)
        a = A(i); b = B(i);
        if any( (pdag(a,:)==-1) & (pdag(:,b)==-1)' );
            pdag(a,b) = -1; pdag(b,a) = 0;
            G(a,b) = 1; G(b,a) = 0;
            %fprintf('rule 2: %d -> %d\n', a, b);
        end
    end
    

    % rule 3
    [A,B] = find(pdag==1); % a-b
    for i=1:length(A)
        a = A(i); b = B(i);
        C = find( (pdag(a,:)==1) & (pdag(:,b)==-1)' );
        % C contains nodes c s.t. a-c->ba
        G2 = setdiag(DAG(C, C), 1);
        if any(G2(:)==0) % there are 2 different non adjacent elements of C
            pdag(a,b) = -1; pdag(b,a) = 0;
            G(a,b) = 1; G(b,a) = 0;
            %fprintf('rule 3: %d -> %d\n', a, b);
        end
    end



end