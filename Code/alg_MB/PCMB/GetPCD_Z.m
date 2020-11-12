
function [PCD,test,time,sepset] = GetPCD_Z(Data,target,alpha,samples,p,maxK)
%
% GetPCD_Z finds the parents, children, and descendants of target node on continuous data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       samples is the number of data samples
%       p is the number of nodes
%       maxK is the maximum size of conditioning set
%
% OUTPUT:
%       PCD is the parents, children, and descendants of the target
%       test is the number of conditional independence tests
%       time is the runtime of the algorithm
%       sepset: the condition sets for each node to make it independt to the target
%
%


if (nargin == 3)
   [samples,p]=size(Data);
   maxK=3;
end

start=tic;


PCD=[];

sepset=cell(1,p);
test=0;

CanPCD = mysetdiff(1:p,target);

last_added=-1;

score=1000000;
    
dep_tmp1 = ones(1,p)*score;
CI=ones(1,p)*score;

dep_tmp_remove = ones(1,p)*score;
CI_remove=ones(1,p)*score;

while( ~isempty(CanPCD) )
    Sep=cell(1,p);
    
    dep=zeros(1,score);
    dep_tmp2 = -score;
    

    tmp_CI=ones(1,p)*score;
    Y=-1;
    
    %----------------GET PC STEP 1-1: get Sep[X]
    
    for i=1:length(CanPCD)
                       
        X = CanPCD(i);
        cutSetSize = 0;
        
        
        while length(PCD) >= cutSetSize &&cutSetSize<=maxK
            
            SS = subsets1(PCD, cutSetSize);   
            
            for si=1:length(SS)
                Z = SS{si};
                
                if last_added~=-1
                    if ~ismember(last_added,Z)
                        continue;
                    end
                end
                
                test=test+1;
                
                [tmp_CI(X),dep(si)]=my_fisherz_test(X,target,Z,Data,samples,alpha);       
                
                if (dep(si) < dep_tmp1(X))          
                    
                    dep_tmp1(X) = dep(si);
                    Sep{X} = Z;
                    CI(X)=tmp_CI(X);
                end
            end

            cutSetSize = cutSetSize + 1;
        end
        
    end
    
    tmp_CanPCD = CanPCD;
    for i=1:length(CanPCD)
                       
        X = CanPCD(i);
        if CI(X)
            sepset{X}=Sep{X};
            tmp_CanPCD = mysetdiff(tmp_CanPCD,X);
        end
    end
    CanPCD = tmp_CanPCD;

    %----------------GET PC STEP 1-1:  END
    
    %----------------GET PC STEP 1-2:  get Y
    
    for i=1:length(CanPCD)
        X = CanPCD(i);
        if dep_tmp1(X)>dep_tmp2&&dep_tmp1(X)~=score
            Y=X;
            dep_tmp2=dep_tmp1(X);
        end
    end
    %----------------GET PC STEP 1-2:  END
        
    %----------------GET PC STEP 1-3:  test CI T with Y 

    if Y~=-1
        if CI(Y)==0||isnan(CI(Y))
            PCD=[PCD Y];
            CanPCD = mysetdiff(CanPCD,Y);
        else
            break;
        end
    else
        break;
    end

    last_added=Y;
    
    %----------------GET PC STEP 1-3:  END
    
    
    % -----------------------------------------------------------
    % remove false positives from PCD
           
    Sep_remove=cell(1,p);
    
    dep_remove=zeros(1,p);   

    tmp_CI_remove=ones(1,p)*score;
    
    for i=1:(length(PCD)-1)

        
                
        X = PCD(i);

        cutSetSize = 0;
        condset = mysetdiff(PCD,X);
               
        while length(condset) >= cutSetSize &&cutSetSize<=maxK
            
            SS = subsets1(condset, cutSetSize);    
            
            for si=1:length(SS)
                Z = SS{si};

                if ~ismember(last_added,Z)
                    continue;
                end
                
                test=test+1;
                [tmp_CI_remove(X),dep_remove(si)]=my_fisherz_test(X,target,Z,Data,samples,alpha);      
                

                if (dep_remove(si) < dep_tmp_remove(X))           
                    
                    dep_tmp_remove(X) = dep_remove(si);
                    Sep_remove{X} = Z;
                    CI_remove(X)=tmp_CI_remove(X);
                end
                
            end
            
            cutSetSize = cutSetSize + 1;
        end
        
    end
    
    tmp_PCD = PCD;
    
    for i=1:(length(PCD)-1)
                        
        X = PCD(i);
        
        if CI_remove(X)
            sepset{X}=Sep_remove{X};
            tmp_PCD = mysetdiff( tmp_PCD,X );

        end
    end
    
    PCD = tmp_PCD;

end

% PCD = sort(PCD);

time=toc(start);





