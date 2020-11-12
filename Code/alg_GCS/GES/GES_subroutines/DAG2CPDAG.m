function Gcp=DAG2CPDAG(G) % transform a DAG to a CPDAG

%% order the edges in G
nodes_order = topo_sort(G); % Perform a topological sort on the nodes of G
% nodes_order(1) is the node which has the highest order
% nodes_order(N) is the node which has the lowest order
edges_order=[];
% edges_order(1,:) is the edge which has the highest order
% edges_order(M,:) is the edge which has the lowest order
M=sum(G(:)); % the number of edges in this DAG
N=size(G,1); % the number of nodes in this DAG

while(size(edges_order,1)<M)
    for ny=N:-1:1
        j=nodes_order(ny);
        inci_all = find(G(:,j)==1)'; % all the edges that incident to j
        if(~isempty(inci_all))
            if(~isempty(edges_order))
                inci = edges_order(find(edges_order(:,2)==j),1)'; % ordered edge that incident to j
                if(~isempty(setdiff(inci_all,inci)))
                    break;
                end
            else
                break;
            end
        end
    end
    for nx=1:N
        i=nodes_order(nx);
        if(~isempty(edges_order))
            if(isempty(find(edges_order(:,2)==j & edges_order(:,1)==i)) & G(i,j)==1)
                break;
            end
        else
            if(G(i,j)==1)
                break;
            end
        end
    end
    edges_order = [edges_order;i,j];
end


%%
sign_edges = zeros(1,M); % 0 means unknown, 1 means compelled, -1 means reversible
while(~isempty(find(sign_edges==0)))
    ss=0;
    for m=M:-1:1 % let x->y be the lowest ordered edge that is labeled "unknown"
        if(sign_edges(m)==0)
            i=edges_order(m,1);
            j=edges_order(m,2);
            break;
        end
    end
    
    idk = find(edges_order(:,2)==i)';
    k= edges_order(idk,1)'; % w->x
    for m=1:length(k)
        if(sign_edges(idk(m))==1)
            if(G(k(m),j)~=1) % if w is not a parent of y
                id=find(edges_order(:,2)==j)'; % label every edge that incident into y with "complled"
                sign_edges(id)=1;
                ss=1;
                break;
            else
                id=find(edges_order(:,1)==k(m) & edges_order(:,2)==j)';% label w->y with "complled"
                sign_edges(id)=1;
            end
        end
    end
    
    if(ss)
        continue;
    end
    
    z=find(G(:,j)==1)';
    if(~isempty(intersect(setdiff(z,i),union(find(G(:,i)==0)',find(G(:,i)==-1)))))
        id = find(edges_order(:,1)==i & edges_order(:,2)==j)';
        sign_edges(id)=1; % label x->y with "compelled"
        
        id1 = find(edges_order(:,2)==j)';
        id2 = intersect(find(sign_edges==0),id1);
        sign_edges(id2)=1; % label all "unknown" edges incident into y with "complled"
    else
        id = find(edges_order(:,1)==i & edges_order(:,2)==j)';
        sign_edges(id)=-1; % label x->y with "reversible"
        
        id1 = find(edges_order(:,2)==j)';
        id2 = intersect(find(sign_edges==0),id1);
        sign_edges(id2)=-1; % label all "unknown" edges incident into y with "reversible"
    end
end

% create CPDAG accoring the labelled edge
Gcp=zeros(N,N);
for m=1:M
    if(sign_edges(m)==1)
        Gcp(edges_order(m,1),edges_order(m,2))=1;
    else
        Gcp(edges_order(m,1),edges_order(m,2))=-1;
        Gcp(edges_order(m,2),edges_order(m,1))=-1;
    end
end










