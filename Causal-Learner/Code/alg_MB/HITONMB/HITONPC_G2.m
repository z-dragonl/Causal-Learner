
function [PC,test,time,sepset] = HITONPC_G2(Data,target,alpha,ns,p,maxK)
%
% HITONPC_G2 finds the parents and children of target node on discrete data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       ns is the size array of each node
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
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end


start=tic;
cpc=[];
PC=[];

dep=zeros(1,p);
sepset=cell(1,p);
test=0;


for i=1:p
    if i~=target
        
        test=test+1;
        [pval, dep(i)]=my_g2_test(i,target,[],Data,ns,alpha);

        if isnan(pval)
            CI=1;
        else
            if pval<=alpha
                CI=0;
            else
                CI=1;
            end
        end
        
        if CI==1
            continue;
        end
        if CI==0
            cpc=[cpc,i];
        end
    end
end

% cpc

[dep_sort,var_index]=sort(dep(cpc),'descend');


for i=1:length(cpc)
    

    Y=cpc(var_index(i));
    PC=[PC Y];
    
    
    pc_length=length(PC);
    pc_tmp=PC;
    
    last_break_flag=0;
    
    for j=pc_length:-1:1
        
        X=PC(j);
        CanPC=mysetdiff(pc_tmp, X);
        
        break_flag=0;
        cutSetSize = 0;
                
        while length(CanPC) >= cutSetSize &&cutSetSize<=maxK
            
            SS = subsets1(CanPC, cutSetSize);   
            
            for si=1:length(SS)
                Z = SS{si};
                
                
                if X~=Y           
                    if isempty(find(Z==Y, 1))
                        continue;
                    end
                end   
                
                
                test=test+1;
                [pval]=my_g2_test(X,target,Z,Data,ns,alpha);      
                if isnan(pval)
                    CI=0;
                else
                   if pval<=alpha
                      CI=0;
                   else
                      CI=1;
                   end
                end
               
                if CI==1
                    pc_tmp=CanPC;
                    sepset{1,X}=Z;
                    break_flag=1;
                    
                    if X==Y
                        last_break_flag=1;
                    end
                    
                    break;
                end
            end
            
            if( break_flag==1 )
                break;
            end
            if( last_break_flag==1 )
                break;
            end
            
            cutSetSize = cutSetSize + 1;
        end

        if( last_break_flag==1 )
            break;
        end
    end
    
    PC=pc_tmp;
    
end

% PC = sort(PC);

time=toc(start);





