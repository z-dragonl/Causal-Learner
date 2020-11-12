function s = check_clique(G, subnode) % check whether node subnode is a clique in G
% here G is a CPDAG
% the definition of clique here: a clique is defined in an undirected graph
% when you ignore the directionality of any directed edges

Gs = G(subnode,subnode); % extract the subgraph
ns = length(subnode);

if(~ns)
    s=1;
else
    [row,col] = find(Gs==1);
    Gs(row,col) = -1;
    Gs(col,row) = -1;
    if((eye(ns)-ones(ns,ns))== Gs) % check whether it is a clique
        s=1;
    else
        s=0;
    end
end