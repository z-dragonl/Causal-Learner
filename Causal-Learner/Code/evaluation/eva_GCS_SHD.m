
function [SHD,reverse,miss,extra,undirected]=eva_GCS_SHD(G,DAG)


diff_G=G-DAG;
[x,y]=find(diff_G~=0);

undirected=0;
reverse=0;
miss=0;
extra=0;
total=0;

for loop_i=1:length(x)
    
    if x(loop_i)==-1
        continue;
    end
    biaozhi=0;
    
    if ~isempty(find(y==(x(loop_i)), 1))
        y_index=find(y==(x(loop_i)));
        for k=1:length(y_index)
            if x(y_index(k))==y(loop_i)
                
                biaozhi=1;
                if G(x(loop_i),y(loop_i))==G(x(y_index(k)),y(y_index(k)))
                    extra=extra+1;
                else
                    reverse=reverse+1;
                    
                end
                x(y_index(k))=-1;
                y(y_index(k))=-1;
                break;
                
            end
        end
    end
    if biaozhi==0
        if G(y(loop_i),x(loop_i))==0
            if G(x(loop_i),y(loop_i))==1
                miss=miss+1;
                
            else
                extra=extra+1;
                
            end
        else
            undirected=undirected+1;
        end
    end
end


total=extra+miss+reverse+undirected;



C=G-DAG;
[x,y]=find(C~=0);

SHD=length(x);

for i=1:length(x)
    if ~isempty(find(y==(x(i))))
        y_index=find(y==(x(i)));
        for k=1:length(y_index)
            if (x(y_index(k)))==(y(i))
                SHD=SHD-0.5;
            end
        end
    end
end

if total~=SHD
    disp('error*****************');
end

