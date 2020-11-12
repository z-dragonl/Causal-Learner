
function [arrhd_F1,arrhd_precision,arrhd_recall,arrhd_distance,true]=eva_GCS_arrhd(G,DAG,nundirected,nreverse,nmiss,nextra,ntotal)



all_true=length(find(G==1));


[x,y]=find(DAG~=0);
all_get=length(x);

for i=1:length(x)
    if ~isempty(find(y==(x(i))))
        y_index=find(y==(x(i)));
        for k=1:length(y_index)
            if (x(y_index(k)))==(y(i))
                all_get=all_get-0.5;
            end
        end
    end
end



inter=all_true-nundirected-nreverse-nmiss;

true=inter;

% all_true
% all_get
% inter

if all_true==0
    if all_get==0
        arrhd_precision=1;
        arrhd_recall=1;
        arrhd_distance=0;
        arrhd_F1=1;
    elseif all_get~=0
        arrhd_precision=0;
        arrhd_recall=0;
        arrhd_distance=sqrt(2);
        arrhd_F1=0;
    end
elseif all_true~=0
    if all_get~=0
        precision_tmp=inter/all_get;
        recall_tmp=inter/all_true;
        distance_tmp=sqrt((1-recall_tmp)*(1-recall_tmp)+(1-precision_tmp)*(1-precision_tmp));
        if (precision_tmp+recall_tmp)==0
            f1_tmp=0;
        else
            f1_tmp=2*precision_tmp*recall_tmp/(precision_tmp+recall_tmp);
        end
        arrhd_precision=precision_tmp;
        arrhd_recall=recall_tmp;
        arrhd_distance=distance_tmp;
        arrhd_F1=f1_tmp;
    else
        arrhd_precision=0;
        arrhd_recall=0;
        arrhd_distance=sqrt(2);
        arrhd_F1=0;
    end
end

