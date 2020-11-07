% V=Delete_validity_test(G, X, Y, H); % do validity test for the operator Delete; V=1 means valid, V=0 mean invalid;
function V = Delete_validity_test(G, i,j, H)
% here G is CPDAG
V=0;

% condition 1
Hj = find(G(j,:)==-1); % neighbors of Xj
Hi = union(find(G(i,:)~=0),find(G(:,i)~=0)'); % adjacent to Xi;
NA = intersect(Hj,Hi); % find the neighbours of Xj and are adjacent to Xi
s1 = check_clique(G,setdiff(NA,H)); % check whether it is a clique

if(s1)
    V=1;
end