function [MB,test] = SRFS_minor_G2(Data,target,alpha, ns, p)
%���㷨Υ����FBED��˼�룬��һ��������ֻ��PC
%����֮�������һ�βü����ü�ʱ���һ���ڵ㲻�ò�
%�ڿռ���������T�����Ľڵ㲻�ٿ��Ǽ���or_rank��ȥ
test=0;
MB=[];
CanMB = mysetdiff([1:p],target);

%��һ��pc����
tag=0;
while(~isempty(CanMB))
    tmp_dep = -10000;
    tmp_pval = 10000;
    tmp_length=length(CanMB);
    
    CanMB_remove=[];
    for i=1:tmp_length
        X = CanMB(i);
        [pval,dep]=my_g2_test(X,target,MB,Data,ns,alpha);
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

%mb�Ĳü�����
tmp_MB = MB;
for i=1:length(MB)-1 %���һ������Ľڵ㲻�ò�
    condset=mysetdiff(MB,MB(i));
    [pval]=my_g2_test(MB(i),target,condset,Data,ns,alpha);
    test = test+1;
    if pval>alpha
        tmp_MB=mysetdiff(tmp_MB,MB(i));
    end
end
MB = tmp_MB;

%�ڶ���Ѱ��T������mb
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
            [pval,dep]=my_g2_test(X,target,MB,Data,ns,alpha);
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

%mb�Ĳü�����
tmp_MB = MB;
for i=1:length(MB)-1 %���һ������Ľڵ㲻�ò�
    condset=mysetdiff(MB,MB(i));
    [pval]=my_g2_test(MB(i),target,condset,Data,ns,alpha);
    test = test+1;
    if pval>alpha
        tmp_MB=mysetdiff(tmp_MB,MB(i));
    end
end

% MB = sort(tmp_MB);
MB = tmp_MB;

end
