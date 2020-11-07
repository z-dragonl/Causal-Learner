function [MB,test,time] = CorrectMB_G2(Data,target,alpha,ns,p,maxK)

start=tic;

MB=[];

test=0;

tmp_CanMB=[];
CanMB=mysetdiff([1:p],target);

while( ~isempty(CanMB) )
    
    %---------------------------------------------------------
    %  add true positives to MB
    
    tmp_dep = -10000;
    tmp_pval = 10000;
    
    for i=1:length(CanMB)
        X = CanMB(i);
        
        test=test+1;
        [pval,dep]=my_g2_test(X,target,MB,Data,ns,alpha);      
        
        if isnan(pval)
            tmp_CanMB = [tmp_CanMB,X];
        else
            if (dep > tmp_dep)           
                
                tmp_dep = dep;
                tmp_pval = pval;
                Y = X;
                
            end
        end
    end
    
    if tmp_pval<=alpha
        MB=[MB Y];
        tmp_CanMB=[tmp_CanMB,Y];
        CanMB= mysetdiff(CanMB,tmp_CanMB);
    else
        
        if tmp_pval>alpha 
            break;
        end
    end
    
    
    % -----------------------------------------------------------
    % remove false positives from MB
    
    tmp_MB = MB;
    for i=1:length(MB)
        condset=mysetdiff(MB,MB(i));
        [pval]=my_g2_test(MB(i),target,condset,Data,ns,alpha);
        test=test+1;
        if pval>alpha
            tmp_MB=mysetdiff(tmp_MB,MB(i));
        end
    end
    
    MB = tmp_MB;
    
end

time=toc(start);





