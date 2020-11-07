% V=Insert_validity_test(G, X, Y, T,1); % do validity test for the operator Insert; V=1 means valid, V=0 mean invalid;
function V = Insert_validity_test2(G, i,j, T)
% here G is CPDAG
V=0;
Tj = find(G(j,:)==-1); % neighbors of Xj
Ti = union(find(G(i,:)~=0),find(G(:,i)~=0)'); % adjacent to Xi;
NA = intersect(Tj,Ti); % find the neighbours of Xj and are adjacent to Xi

% condition 2: every semi-directed path from Xj to Xi contains a node in union(NA,T)
% Note: EVERY!!
s2=Insert_vC2_new(G,j,i,union(NA,T));
if(s2)
    V=1;
end


