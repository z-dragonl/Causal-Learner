%calculate the changed score after the insert operator: i->j
function [chscore,desc,record_local_score] = Insert_changed_score(Data,G,i,j,T,record_local_score,score_func,parameters)

Tj = find(G(j,:)==-1); % neighbors of Xj
Ti = union(find(G(i,:)~=0),find(G(:,i)~=0)'); % adjacent to Xi;
NA = intersect(Tj,Ti); % find the neighbours of Xj and are adjacent to Xi
Paj = find(G(:,j)==1)'; % find the parents of Xj
% the function local_score() calculates the local score
tmp1=union(NA,T);
tmp2=union(tmp1,Paj);
tmp3=union(tmp2,i);

% before you calculate the local score, firstly you search in the
% "record_local_score", to see whether you have calculated it before
r=length(record_local_score{j});
s1=0;s2=0;

for r0 = 1:r
    if(isequal(record_local_score{j}{r0}(1:end-1),tmp3)) 
        score1 = record_local_score{j}{r0}(end);
        s1=1;
    end
    
    if(isequal(record_local_score{j}{r0}(1:end-1),tmp2)) % notice the differnece between 0*0 empty matrix and 1*0 empty matrix
        score2 = record_local_score{j}{r0}(end);
        s2=1;
    else
        if(isequal(record_local_score{j}{r0}(1:end-1),0) && isempty(tmp2))
            score2 = record_local_score{j}{r0}(end);
            s2=1;
        end
    end
    
    if(s1 & s2)
        break;
    end
end

if(~s1)
    score1 = feval(score_func,Data,j,tmp3,parameters);
    record_local_score{j}{r+1}=[tmp3,score1];
end
if(~s2)
    score2 = feval(score_func,Data,j,tmp2,parameters);
    r=length(record_local_score{j});
    if(~isempty(tmp2))
        record_local_score{j}{r+1}=[tmp2,score2];
    else
        record_local_score{j}{r+1}=[0,score2];
    end
end
chscore = score1-score2;
desc{1}=i;
desc{2}=j;
desc{3}=T;








