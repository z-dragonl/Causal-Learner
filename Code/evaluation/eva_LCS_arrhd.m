
function [arrhd_F1,arrhd_precision,arrhd_recall,arrhd_distance]=eva_LCS_arrhd(PC,P,C,trueP,trueC,truePC)


if isempty(truePC)
    if isempty(PC)
        arrhd_precision=1;
        arrhd_recall=1;
        arrhd_distance=0;
        arrhd_F1=1;
    elseif ~isempty(PC)
        arrhd_precision=0;
        arrhd_recall=0;
        arrhd_distance=sqrt(2);
        arrhd_F1=0;
    end
elseif ~isempty(truePC)
    if ~isempty(PC)
        precision_tmp=(length(intersect(P,trueP))+length(intersect(C,trueC)))/length(PC);
        recall_tmp=(length(intersect(P,trueP))+length(intersect(C,trueC)))/length(truePC);
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

