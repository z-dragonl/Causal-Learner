function s = check2(G,Nx,Ax) % check for every nodes in Nx, it is adjacent to all the nodes in Ax except itself

s=1;
for i=1:length(Nx)
    j = setdiff(Ax,Nx(i));
    if(~isempty(intersect(find(G(Nx(i),j)==0),find(G(j,Nx(i))'==0))))
        s=0;
        break;
    end
end

