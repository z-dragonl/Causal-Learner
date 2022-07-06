function [MB,test] = SRFS_minor_Z(Data,target,alpha, ns, p)
%该算法违背了FBED的思想，第一轮理论上只找PC
%两轮之间多添加一次裁剪，裁剪时最后一个节点不用测
%在空集条件下与T独立的节点不再考虑加入or_rank中去
test=0;
MB=[];
[samples,p]=size(Data);
CanMB = mysetdiff([1:p],target);

%第一轮pc发现
tag=0;
while(~isempty(CanMB))
    tmp_dep = -10000;
    tmp_pval = 10000;
    tmp_length=length(CanMB);
    
    CanMB_remove=[];
    for i=1:tmp_length
        X = CanMB(i);
        [pval,dep]=my_fisherz_test(X,target,MB,Data,samples,alpha);
        test = test+1;
        if ~isnan(pval)
            if (dep > tmp_dep)
                tmp_dep = dep;
                tmp_pval = pval;
                Y = X;
            end
            if pval>alpha
                CanMB_remove=[CanMB_remove X];
            end
        else
            CanMB_remove=[CanMB_remove X];
        end
    end
    CanMB=mysetdiff(CanMB,CanMB_remove);
    
    if tmp_pval<=alpha
        MB=[MB Y];
        tag=1;
        CanMB=mysetdiff(CanMB,Y);
    end
end

%mb的裁剪操作
tmp_MB = MB;
for i=1:length(MB)-1 %最后一个加入的节点不用测
    condset=mysetdiff(MB,MB(i));
    [pval]=my_fisherz_test(MB(i),target,condset,Data,samples,alpha);
    test = test+1;
    if pval>alpha
        tmp_MB=mysetdiff(tmp_MB,MB(i));
    end
end
MB = tmp_MB;

%第二轮寻找T的所有mb
if tag==1
    CanMB = mysetdiff([1:p],target);
    CanMB = mysetdiff(CanMB,MB);
    while(~isempty(CanMB))
        tmp_dep = -10000;
        tmp_pval = 10000;
        tmp_length=length(CanMB);
        
        CanMB_remove=[];
        for i=1:tmp_length
            X = CanMB(i);
            [pval,dep]=my_fisherz_test(X,target,MB,Data,samples,alpha);
            test=test+1;
            if ~isnan(pval)
                if (dep > tmp_dep)
                    tmp_dep = dep;
                    tmp_pval = pval;
                    Y = X;
                end
                if pval>alpha
                    CanMB_remove=[CanMB_remove X];
                end
            else
                CanMB_remove=[CanMB_remove X];
            end
        end
        CanMB=mysetdiff(CanMB,CanMB_remove);
        
        if tmp_pval<=alpha
            MB=[MB Y];
            CanMB=mysetdiff(CanMB,Y);
        end
    end
end

%mb的裁剪操作
tmp_MB = MB;
for i=1:length(MB)-1 %最后一个加入的节点不用测
    condset=mysetdiff(MB,MB(i));
    [pval]=my_fisherz_test(MB(i),target,condset,Data,samples,alpha);
    test = test+1;
    if pval>alpha
        tmp_MB=mysetdiff(tmp_MB,MB(i));
    end
end

% MB = sort(tmp_MB);
MB = tmp_MB;

end

