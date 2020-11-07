function [MB,test,time] = IAMB_Forward_G2(Data,target,alpha,ns,p,maxK)

% Forward part of IAMB_G2

start=tic;


MB=[];

test=0;

tmp_CanMB=[];
CanMB = mysetdiff([1:p],target);

%---------------------------------------------------------
%  add true positives to MB

while( ~isempty(CanMB) )
    
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
    else
        if tmp_pval>alpha 
            break;
        end
    end
    
    CanMB= mysetdiff(CanMB,tmp_CanMB);
end


time=toc(start);





