function nodes_order = topo_sort(G)  % Perform a topological sort on the nodes of G
% Note here G is a DAG
%topological sort refers to total ordering of the nodes where if Xi is
%ancestor of Xj, then Xi must precede Xj in the ordering

N=size(G,1);
nodes_order=[]; % nodes_order(1) is the node with the highest order, nodes_order(N) is the node with the lowest order
sign=zeros(1,N); % sign(i)=0 means the node has not been ordered
while(length(nodes_order)<N)
    for i=1:N
        if(~sign(i))
            index=find(G(:,i)==1)'; % find the parents of i
            index0=find(sign==0);
            if(isempty(intersect(index,index0)))
                nodes_order = [nodes_order,i];
                sign(i)=1;
                break;
            end
        end
    end
end



