
function [adj_F1,adj_precision,adj_recall,adj_distance]=eva_MB_adj(MB,trueMB)



if isempty(trueMB)
    if isempty(MB)
        adj_precision=1;
        adj_recall=1;
        adj_distance=0;
        adj_F1=1;
    elseif ~isempty(MB)
        adj_precision=0;
        adj_recall=0;
        adj_distance=sqrt(2);
        adj_F1=0;
    end
elseif ~isempty(trueMB)
    if ~isempty(MB)
        precision_tmp=length(intersect(MB,trueMB))/length(MB);
        recall_tmp=length(intersect(MB,trueMB))/length(trueMB);
        distance_tmp=sqrt((1-recall_tmp)*(1-recall_tmp)+(1-precision_tmp)*(1-precision_tmp));
        if (precision_tmp+recall_tmp)==0
            f1_tmp=0;
        else
            f1_tmp=2*precision_tmp*recall_tmp/(precision_tmp+recall_tmp);
        end
        adj_precision=precision_tmp;
        adj_recall=recall_tmp;
        adj_distance=distance_tmp;
        adj_F1=f1_tmp;
    else
        adj_precision=0;
        adj_recall=0;
        adj_distance=sqrt(2);
        adj_F1=0;
    end
end

