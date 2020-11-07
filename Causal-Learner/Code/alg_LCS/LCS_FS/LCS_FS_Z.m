
function [P,C,PC,UN,test,time] = LCS_FS_Z(Data,target,alpha,samples,p,maxK)
%
% LCS_FS_Z finds and distinguishes the parents and children of target node on continuous data
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
%       P is the parents of the target
%       C is the children of the target
%       PC is the union of the parents and children of the target
%       UN is the nodes in PC but cannot distinguish whether they are parents or children
%       test is the number of conditional independence tests
%       time is the runtime of the algorithm
%
%


if (nargin == 3)
   [samples,p]=size(Data);
   maxK=3;
end
                                   

start=tic;

test=0;

all_sepset=cell(1,p);
DAG=zeros(p,p);
pdag=zeros(p,p);
G=zeros(p,p);

PC_calculated=ones(1,p);

all_PC=cell(1,p);

Q=target;
Tmp=[];

while length(Tmp)<=p && ~isempty(Q)

    
    A=Q(1);    Q(1)=[];
    
    if ismember(A,Tmp)
        continue;
    else
        Tmp=[Tmp A];
    end
    
    
    if ~ismember(1,pdag(A,:))&&A~=target
        continue;
    end
    
    if PC_calculated(A)
        all_PC{A}=FCBF_PC_Z(Data,A,alpha,samples,p,maxK);
        PC_calculated(A)=0;
    end
    
    for i=1:length(all_PC{A})
        
        B=all_PC{A}(i);   Q=[Q B];
        
        % create undirected edges between A and B£º A - B 
        
        DAG(A,B)=1;  DAG(B,A)=1;
        
        if pdag(A,B)==0 && pdag(B,A)==0
            pdag(A,B)=1;  pdag(B,A)=1;
            G(A,B)=1;  G(B,A)=1;
        end
        
        if PC_calculated(B)
            [all_PC{B},ntest_PC]=FCBF_PC_Z(Data,B,alpha,samples,p,maxK);
            PC_calculated(B)=0;
            
            test=test+ntest_PC;
        end
        
        for j=1:length(all_PC{B})
            
            C=all_PC{B}(j);
            
            % create undirected edges between B and C£º B - C 
            
            DAG(C,B)=1;  DAG(B,C)=1;      
            
            if pdag(C,B)==0 && pdag(B,C)==0
                pdag(C,B)=1;  pdag(B,C)=1;
                G(C,B)=1;  G(B,C)=1;
            end
            
            if ismember(C,all_PC{A})||C==A
                continue;
            end
            
            CanPC=all_PC{A};
            cutSetSize = 0;
            
            break_flag=0;
            
            while length(CanPC) >= cutSetSize &&cutSetSize<=maxK
                
                SS = subsets1(CanPC, cutSetSize);   
                for si=1:length(SS)
                    Z = SS{si};
                    
                    test=test+1;
                    CI=my_fisherz_test(C,A,Z,Data,samples,alpha);
                    
                    if CI
                        
                        all_sepset{A}{C}=Z;
                        
                        if ~ismember(B,Z)
                            test=test+1;
                            CI=my_fisherz_test(C,A,myunion(Z,B),Data,samples,alpha);
                            
                            if isnan(CI)||CI==0
%                                 fprintf('V-structure: %d->%d<-%d\n',A,B,C);
                                
                                pdag(A,B) = -1; pdag(C,B) = -1; pdag(B,A) = 0; pdag(B,C) = 0;
                                G(A,B) = 1;     G(C,B) = 1;     G(B,A) = 0;    G(B,C) = 0;

                            end
                        end
                        
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
    end
    
    % Meek rules
    [DAG,pdag,G]=meeks(DAG,pdag,G,p);
    
    if ~ismember(1,pdag(target,:))&&~ismember(1,pdag(:,target))
        % P and C have been distinguished
        break;
    end
    
end

P=find(pdag(:,target)==-1);
C=find(pdag(target,:)==-1);
UN=find(pdag(target,:)==1);

PC=myunion(myunion(P,C),UN);

time=toc(start);


