
function [MB,test,time] = BAMB_G2(Data,target,alpha,ns,p,maxK)
%
% BAMB_G2 finds the Markov blanket of target node on discrete data
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
%       MB is the Markov blanket of the target
%       test is the number of conditional independence tests
%       time is the runtime of the algorithm
%
%


if (nargin == 3)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end


start=tic;

test=0;

cpc=[];
pc=[];


dep=zeros(1,p);
sepset=cell(1,p);

spouse = cell(1,p);

dep_sp=zeros(p,p);

% Remove the nodes independent of the target node conditioning on
% empty set, and sort the dependent nodes

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

[~,var_index]=sort(dep(cpc),'descend');

CI_dependent=[];
for i=1:length(cpc)
    CI_dependent=[CI_dependent cpc(var_index(i))];
end


for i=1:length(cpc)
    pc=[pc cpc(var_index(i))];
    pc_length=length(pc);
    pc_tmp=pc;
    last_break_flag=0;
    
    % Step 1: Find the candidate set of PC and candidate set of spouses
    
    for j=pc_length:-1:1
        Y = pc(j);
        
        CanPC=mysetdiff(pc_tmp, Y);
        
        cutSetSize = 1;            
        other_break_flag=0;
        while length(CanPC) >= cutSetSize && cutSetSize<=maxK
            SS = subsets1(CanPC, cutSetSize);   
            for si=1:length(SS)
                Z = SS{si};
                
                if cpc(var_index(i))~=Y
                    if isempty(find(Z==cpc(var_index(i)), 1))
                        continue;
                    end
                end
                
                
                
                test=test+1;
                [pval]=my_g2_test(Y,target,Z,Data,ns,alpha);       
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
                    sepset{Y}=Z;
                    
                    spouse{1,Y}=[];        
                    
                    if cpc(var_index(i))==Y
                        last_break_flag=1;
                    end
                    other_break_flag=1;  
                    break;
                end
            end
            
            if( last_break_flag==1 )
                break;
            end
            if( other_break_flag==1 )
                break;
            end
            cutSetSize = cutSetSize + 1;
        end
        
        if( last_break_flag==1 )
            break;
        end
        
        
        if cpc(var_index(i))==Y
            
            
            NonTPC = mysetdiff(1:p,myunion(pc_tmp,target));
            
            for k=1:length(NonTPC)
                X = NonTPC(k);
                
                S =  myunion(sepset{1,X}, Y);      
                test = test +1;
                [pval,dep_sp(Y,X)]=my_g2_test(target,X,S,Data,ns,alpha);      
                if isnan(pval)
                    CI=1;
                else
                    if pval<=alpha
                        CI=0;
                    else
                        CI=1;
                    end
                end
                
                if CI==0
                    spouse{1,Y} = myunion ( spouse{1,Y} ,X );
                end
            end
            
            
            % Step 2: Remove false positives from the candidate set of spouses
            % Step 2-1
            
            SP=[];
            [~,SP_index]=sort(dep_sp(Y,spouse{1,Y}),'descend');
            for f=1:length(spouse{1,Y})
                SP=[SP spouse{1,Y}(SP_index(f))];
                SP_length=length(SP);
                SP_tmp=SP;
                last_SP_break_flag1=0;
                for c=SP_length:-1:1
                    X = SP(c);
                    CanSP=mysetdiff(SP_tmp, X);
                    cutSetSize = 1;           
                    SS = subsets1(CanSP, cutSetSize);   
                    for si=1:length(SS)
                        Z = SS{si};
                        condset=myunion(Z,Y);
                        
                        if spouse{1,Y}(SP_index(f))~=X
                            if isempty(find(Z==spouse{1,Y}(SP_index(f)), 1))
                                continue;
                            end
                        end
                        
                        test=test+1;
                        [pval]=my_g2_test(X,target,condset,Data,ns,alpha);
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
                            SP_tmp=CanSP;
                            
                            if spouse{1,Y}(SP_index(f))==X
                                last_SP_break_flag1=1;
                            end
                            break;
                        end
                    end
                    if( last_SP_break_flag1==1 )
                        break;
                    end
                end
                SP=SP_tmp;
            end
            spouse{1,Y} = SP;
            
            
            % Step 2-2
            
            SP=[];
            for f=1:length(spouse{1,Y})
                SP=[SP spouse{1,Y}(f)];
                SP_length=length(SP);
                SP_tmp=SP;
                last_SP_break_flag2=0;
                for c=SP_length:-1:1
                    X = SP(c);
                    CanSP=mysetdiff(SP_tmp, X);
                    cutSetSize = 0;           
                    
                    other_SP_break_flag2=0;
                    while length(CanSP) >= cutSetSize && cutSetSize<=maxK
                        SS = subsets1(CanSP, cutSetSize);   
                        for si=1:length(SS)
                            Z = SS{si};
                            
                            condset=Z;
                            
                            if spouse{1,Y}(f)~=X
                                if isempty(find(Z==spouse{1,Y}(f), 1))
                                    continue;
                                end
                            end
                            
                            test=test+1;
                            [pval]=my_g2_test(X,Y,condset,Data,ns,alpha);
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
                                
                                SP_tmp=mysetdiff(SP_tmp, X);
                                
                                if spouse{1,Y}(f)==X
                                    last_SP_break_flag2=1;
                                end
                                other_SP_break_flag2=1;
                                break;
                            end
                        end
                        if last_SP_break_flag2||other_SP_break_flag2
                            break;
                        end
                        cutSetSize = cutSetSize + 1;
                    end
                    if last_SP_break_flag2
                        break;
                    end
                end
                SP=SP_tmp;
            end
            spouse{1,Y} = SP;
            
            
            % Step 2-3
            
            SP=[];
            for f=1:length(spouse{1,Y})
                SP=[SP spouse{1,Y}(f)];
                SP_length=length(SP);
                SP_tmp=SP;
                last_SP_break_flag3=0;
                for c=SP_length:-1:1
                    X = SP(c);
                    CanSP=mysetdiff(myunion(pc_tmp,SP_tmp), X);
                    cutSetSize = 1;          
                    
                    other_SP_break_flag3=0;
                    while length(CanSP) >= cutSetSize && cutSetSize<=maxK
                        SS = subsets1(CanSP, cutSetSize);  
                        for si=1:length(SS)
                            Z = SS{si};
                            
                            condset=unique([Z,Y]);
                            
                            if (length(intersect(pc,condset))==1)
                                continue;
                            end
                            
                            
                            if spouse{1,Y}(f)~=X
                                if isempty(find(Z==spouse{1,Y}(f), 1))
                                    continue;
                                end
                            end
                            
                            test=test+1;
                            [pval]=my_g2_test(X,target,condset,Data,ns,alpha);
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
                                
                                SP_tmp=mysetdiff(SP_tmp, X);
                                
                                if spouse{1,Y}(f)==X
                                    last_SP_break_flag3=1;
                                end
                                other_SP_break_flag3=1;
                                break;
                            end
                        end
                        if last_SP_break_flag3||other_SP_break_flag3
                            break;
                        end
                        cutSetSize = cutSetSize + 1;
                    end
                    if last_SP_break_flag3
                        break;
                    end
                end
                SP=SP_tmp;
            end
            spouse{1,Y} = SP;
            
            
            % Step 3: Remove false positives from the candidate set of PC
            
            M=pc_tmp;
            
            pc_length=length(M);
            pc_Tmp=M;
            last_PC_break_flag=0;
            
            for long=pc_length:-1:1
                Y = M(long);
                CanPC=mysetdiff(pc_Tmp, Y);
                
                cutSetSize = 1;            
                other_PC_break_flag=0;
                while length(CanPC) >= cutSetSize && cutSetSize<=maxK
                    SS = subsets1(CanPC, cutSetSize);   
                    for si=1:length(SS)
                        Z = SS{si};
                        
                        if cpc(var_index(i))~=Y
                            if isempty(find(Z==cpc(var_index(i)), 1))
                                continue;
                            end
                        end
                        
                        
                        spouse_test=spouse{1,Z};
                        
                        if isempty(spouse_test)
                            continue;
                        end
                        
                        TestSet = myunion(Z,spouse_test);
                        
                        test=test+1;
                        [pval]=my_g2_test(Y,target,TestSet,Data,ns,alpha);      
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
                            M=mysetdiff( M,Y );
                            spouse{1,Y} = [];
                            if cpc(var_index(i))==Y
                                last_PC_break_flag=1; 
                            end
                            other_PC_break_flag=1;     
                            break;
                        end
                    end
                    if other_PC_break_flag||last_PC_break_flag
                        break;
                    end
                    cutSetSize = cutSetSize + 1;
                end
                if last_PC_break_flag
                    break;
                end
            end

            pc_tmp=M;

        end
    end
    
    pc=pc_tmp;
end


MB=myunion(pc,cell2mat(spouse));

time=toc(start);


