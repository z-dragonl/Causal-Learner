function Gd = PDAG2DAG(G) % transform a PDAG to DAG

% first create a DAG that contains all the directed edges in PDAG
Gd=G;
for i=1:size(G,1)
    for j=1:size(G,2)
        Gd(i,j)=max(Gd(i,j),0);
    end
end

Gp=G;
inde=zeros(1,size(Gp,1)); % index whether the ith node has been removed. 1:removed; 0: not
% while(~isempty(find(inde==0)))
%     for i=1:size(Gp,1)
%         if(~inde(i))
%             sign=0;
%             if(isempty(intersect(find(Gp(i,:)==1), find(inde==0)))) % Xi has no out-going edges
%                 sign=sign+1;
%                 Nx=intersect(find(Gp(i,:)==-1), find(inde==0)); % find the neighbors of Xi in P
%                 Pax = intersect(find(Gp(:,i)==1)', find(inde==0)); % find the parents of Xi in P
%                 if(~isempty(Nx))
%                     if(check_clique(Gp,union(Nx,Pax))) % check whether it is a clique   Here May have problems!!!!
%                         sign=sign+1;
%                     end
%                 else
%                     sign=sign+1;
%                 end
%             end
%             if(sign==2)
%                 Gd(intersect(find(Gp(i,:)==-1),find(inde==0)),i) = 1; % for each undirected edge Y-X in PDAG, insert a directed edge Y->X in G
%                 Gd(i,intersect(find(Gp(i,:)==-1),find(inde==0))) = 0;
%                 inde(i)=1;
%                 %             Gp=[Gp(:,1:i-1),Gp(:,i+1:end)]; % remove Xi. can't remove directly which will lead to
%                 %             Gp=[Gp(1:i-1,:);Gp(i+1:end,:)];
%             end
%         end
%     end
% end

while(~isempty(find(inde==0)))
    for i=1:size(Gp,1)
        if(~inde(i))
            sign=0;
            if(isempty(intersect(find(Gp(i,:)==1), find(inde==0)))) % Xi has no out-going edges
                sign=sign+1;
                Nx=intersect(find(Gp(i,:)==-1), find(inde==0)); % find the neighbors of Xi in P
                Ax = intersect(union(find(Gp(:,i)==1)', find(Gp(i,:)==1)), find(inde==0)); % find the adjacent of Xi in P
                Ax = union(Ax,Nx);
                if(~isempty(Nx))
                    if(check2(Gp,Nx,Ax)) %%% according to the original paper
                        sign=sign+1;
                    end
                else
                    sign=sign+1;
                end
            end
            if(sign==2)
                Gd(intersect(find(Gp(i,:)==-1),find(inde==0)),i) = 1; % for each undirected edge Y-X in PDAG, insert a directed edge Y->X in G
                Gd(i,intersect(find(Gp(i,:)==-1),find(inde==0))) = 0;
                inde(i)=1;
                %             Gp=[Gp(:,1:i-1),Gp(:,i+1:end)]; % remove Xi. can't remove directly which will lead to
                %             Gp=[Gp(1:i-1,:);Gp(i+1:end,:)];
            end
        end
    end
end



