function [MB,ntest,time,index] = SRFS_main_G2(Data,target,alpha, ns, p)
%该算法违背了FBED的思想，第一轮理论上只找PC
%两轮之间多添加一次裁剪，裁剪时最后一个节点不用测
%在空集条件下与T独立的节点不再考虑加入or_rank中去
start=tic;

p_rank=zeros(1,p);%记录每个节点与T之间的依赖度
dele_MB=[];%记录裁剪阶段被删除的节点

MB=[];
ntest=0;
CanMB = mysetdiff([1:p],target);

%第一轮pc发现
CanSP=[];
tag=0;
while(~isempty(CanMB))
    tmp_dep = -10000;
    tmp_pval = 10000;
    tmp_length=length(CanMB);
    
    CanMB_remove=[];
    for i=1:tmp_length
        X = CanMB(i);
        ntest=ntest+1;
        [pval,dep]=my_g2_test(X,target,MB,Data,ns,alpha);
        
        if ~isnan(pval)
            if (dep > tmp_dep)
                tmp_dep = dep;
                tmp_pval = pval;
                Y = X;
            end
            if pval>alpha
                CanMB_remove=[CanMB_remove X];
                if isempty(MB)
                    CanSP=[CanSP X];
                end
            end
            p_rank(X)=pval;
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
    [pval]=my_g2_test(MB(i),target,condset,Data,ns,alpha);
    ntest=ntest+1;
    if pval>alpha
        tmp_MB=mysetdiff(tmp_MB,MB(i));
        dele_MB=[dele_MB MB(i)];
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
            ntest=ntest+1;
            [pval,dep]=my_g2_test(X,target,MB,Data,ns,alpha);
            
            if ~isnan(pval)
                if (dep > tmp_dep)
                    tmp_dep = dep;
                    tmp_pval = pval;
                    Y = X;
                end
                if pval>alpha
                    CanMB_remove=[CanMB_remove X];
                end
                p_rank(X)=pval;
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
    [pval]=my_g2_test(MB(i),target,condset,Data,ns,alpha);
    ntest=ntest+1;
    if pval>alpha
        tmp_MB=mysetdiff(tmp_MB,MB(i));
        dele_MB=[dele_MB MB(i)];
    end
end

% MB = sort(tmp_MB);
MB = tmp_MB;

[p_value,index]=sort(p_rank,'ascend');

index=mysetdiff(index,target);
index=mysetdiff(index,MB);
index=mysetdiff(index,dele_MB);
index=mysetdiff(index,CanSP);

time=toc(start);

end

