function G = Delete(G,i,j,H) % Delete operator

G(i,j)=0; % delete the edge between Xi and Xj
G(j,i)=0;
for k=1:length(H) % directing the previous undirected edge
    G(j,H(k))=1;
    G(H(k),j)=0;
    G(i,H(k))=1;
    G(H(k),i)=0;
end