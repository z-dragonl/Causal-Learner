function [MB,test,time] = IAMB_Forward_Z(Data,target,alpha,samples,p,maxK)

% Forward part of IAMB_Z

start=tic;


MB=[];

test=0;

tmp_CanMB=[];
CanMB = mysetdiff([1:p],target);

%---------------------------------------------------------
%  add true positives to MB

while( ~isempty(CanMB) )
    
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
    else
        if tmp_CI 
            break;
        end
    end
    
    CanMB= mysetdiff(CanMB,tmp_CanMB);
end


time=toc(start);





