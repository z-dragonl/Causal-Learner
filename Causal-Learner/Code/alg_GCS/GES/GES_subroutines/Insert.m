function G = Insert(G,i,j,T) % Insert operator

G(i,j)=1; % insert the directed edge Xi->Xj
for k=1:length(T) % directing the previous undirected edge between T and Xj as T->Xj
    G(T(k),j)=1;
    G(j,T(k))=0;
end