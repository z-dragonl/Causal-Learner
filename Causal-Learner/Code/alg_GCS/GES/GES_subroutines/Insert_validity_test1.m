% V=Insert_validity_test1(G, X, Y, T,1); % do validity test for the operator Insert; V=1 means valid, V=0 mean invalid;
function V = Insert_validity_test1(G, i,j, T)
% here G is CPDAG
V=0;

% condition 1
Tj = find(G(j,:)==-1); % neighbors of Xj
Ti = union(find(G(i,:)~=0),find(G(:,i)~=0)'); % adjacent to Xi;
NA = intersect(Tj,Ti); % find the neighbours of Xj and are adjacent to Xi
V = check_clique(G,union(NA,T)); % check whether it is a clique

