
function [PC,test,time,sepset] = MMPC_Z(Data,target,alpha,samples,p,maxK)
%
% MMPC_Z finds the parents and children of target node on continuous data
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
%       PC is the parents and children of the target
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

PC=[];

sepset=cell(1,p);
test=0;

U = mysetdiff(1:p,target);

score=1000000;

last_added=-1;

CI=ones(1,p)*score;
dep_tmp1 = ones(1,p)*score;

while( ~isempty(U) )
    Sep=cell(1,p);
    

    dep_tmp2 = (-1)*score;
    
    tmp_U = U;
    

    Y=-1;
    
    tmp_CI=ones(1,p)*score;
    dep=zeros(1,score);

    %----------------GET PC STEP 1-1: get Sep[X]  
    
    
    for i=1:length(U)
                       
        X = U(i);
        cutSetSize = 0;
        
        break_flag=0;
        
        while length(PC) >= cutSetSize&&cutSetSize<=maxK
            
            SS = subsets1(PC, cutSetSize);   
            
            for si=1:length(SS)
                Z = SS{si};
                
                if last_added~=-1
                    if ~ismember(last_added,Z)
                        continue;
                    end
                end
                
                test=test+1;
                
                [tmp_CI(X),dep(si)]=my_fisherz_test(X,target,Z,Data,samples,alpha);       
                
                if isnan(tmp_CI(X))
                    tmp_CI(X)=0;
                end
                
                if tmp_CI(X)
                    sepset{X}=Z;
                    tmp_U = mysetdiff(tmp_U,X);
                    
                    break_flag=1;
                    break;
                end
                
                if tmp_CI(X) == 0
                    if (dep(si) < dep_tmp1(X))           
                        
                        dep_tmp1(X) = dep(si);
                        Sep{X} = Z;
                        CI(X)=tmp_CI(X);
                    end
                end
            end
           
            if break_flag
                break;
            end
            
            cutSetSize = cutSetSize + 1;
        end
        
    end
    U =  tmp_U;
    
    %----------------GET PC STEP 1-1:  END
    
    
    %----------------GET PC STEP 1-2:  get Y
    
    for i=1:length(U)
        X = U(i);
        if dep_tmp1(X)>dep_tmp2&&dep_tmp1(X)~=score
            Y=X;
            dep_tmp2=dep_tmp1(X);
        end
    end
    
    %----------------GET PC STEP 1-2:  END
        
    
    %----------------GET PC STEP 1-3:  test CI T with Y 

    if Y~=-1&&CI(Y)==0
        PC=[PC Y];
        U = mysetdiff(U,Y);
    else
        break;
    end
    
    %----------------GET PC STEP 1-3:  END
    
    last_added=Y;
end

% ------------------------------
% remove false positives from PC 

tmp_PC = PC;
for i=(length(PC)-1):-1:1
    break_flag=0;
        
    X = PC(i);

    CanPC=mysetdiff(tmp_PC, X);
        
    cutSetSize = 0;
    
    while length(CanPC) >= cutSetSize&&cutSetSize<=maxK
        
        SS = subsets1(CanPC, cutSetSize);   
        
        for si=1:length(SS)
            Z = SS{si};
                        
            test=test+1;
            
            [CI]=my_fisherz_test(X,target,Z,Data,samples,alpha);   
            
            
            if isnan(CI)
                CI=0;
            end
            
            if CI==1
                tmp_PC = mysetdiff(tmp_PC,X);
                sepset{X}=Z;
                
                break_flag=1;
                break;
            end
        end
        
        if break_flag
            break;
        end
        cutSetSize = cutSetSize + 1;
    end
    
end

PC=tmp_PC;

% PC = sort(tmp_PC);

time=toc(start);





