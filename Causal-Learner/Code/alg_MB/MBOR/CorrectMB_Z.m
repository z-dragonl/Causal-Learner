function [MB,test,time] = CorrectMB_Z(Data,target,alpha,samples,p,maxK)

start=tic;

MB=[];

test=0;

tmp_CanMB=[];
CanMB=mysetdiff([1:p],target);

while( ~isempty(CanMB) )
    
    %---------------------------------------------------------
    %  add true positives to MB
    
    tmp_dep = -10000;
    tmp_CI = 10000;
    
    for i=1:length(CanMB)
        X = CanMB(i);
        
        test=test+1;
        [CI,dep]=my_fisherz_test(X,target,MB,Data,samples,alpha);      
        
        if isnan(CI)
            tmp_CanMB = [tmp_CanMB,X];
        else
            if (dep > tmp_dep)             
                
                tmp_dep = dep;
                tmp_CI = CI;
                Y = X;
                
            end
        end
    end
    
    if tmp_CI==0
        MB=[MB Y];
        tmp_CanMB=[tmp_CanMB,Y];
        CanMB= mysetdiff(CanMB,tmp_CanMB);
    else
        
        if tmp_CI
            break;
        end
    end
    
    
    % -----------------------------------------------------------
    % remove false positives from MB
    
    tmp_MB = MB;
    for i=1:length(MB)
        condset=mysetdiff(MB,MB(i));
        [CI]=my_fisherz_test(MB(i),target,condset,Data,samples,alpha);
        test=test+1;
        if CI
            tmp_MB=mysetdiff(tmp_MB,MB(i));
        end
    end
    
    MB = tmp_MB;
    
end

time=toc(start);





