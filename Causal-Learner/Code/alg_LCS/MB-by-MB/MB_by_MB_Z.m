
function [P,C,PC,UN,test,time] = MB_by_MB_Z(Data,target,alpha,samples,p,maxK)
%
% MB_by_MB_Z finds and distinguishes the parents and children of target node on continuous data
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


DAG=zeros(p,p);
pdag=zeros(p,p);
G=zeros(p,p);

MB_calculated=ones(1,p);

all_PC=cell(1,p);
all_MB=cell(1,p);
all_can_spouse=cell(1,p);
all_sepset=cell(p,p);

Q=target;
Tmp=[];

num_calculated=0;

while length(Tmp)<=p && ~isempty(Q)
    
    
    A=Q(1);    Q(1)=[];
    
    if ismember(A,Tmp)
        continue;
    else
        Tmp=[Tmp A];
    end
    
    
    % Get MB(A)
    if MB_calculated(A)
        [all_MB{A},ntest1]=IAMB_Forward_Z (Data,A,alpha,samples,p,maxK);
        test=test+ntest1;
        MB_calculated(A)=0; 
    end

    
    all_PC{A}=all_MB{A};
    

    for i=1:length(all_MB{A})
        
        B=all_MB{A}(i);   Q=[Q B];
        
        DAG(A,B)=1;  DAG(B,A)=1;       
        
        if pdag(A,B)==0 && pdag(B,A)==0
            pdag(A,B)=1;  pdag(B,A)=1;
            G(A,B)=1;     G(B,A)=1;
        end
        
        % Find PC
        cutSetSize = 0;
        break_flag=0;
        CanPC=mysetdiff(all_MB{A},B);
        while length(CanPC) >= cutSetSize &&cutSetSize<=maxK
            SS = subsets1(CanPC, cutSetSize);  
            for si=1:length(SS)
                Z = SS{si};
                test=test+1;
                CI=my_fisherz_test(B,A,Z,Data,samples,alpha);
                
                if CI
                    all_sepset{A}{B}=Z;  all_sepset{B}{A}=Z;
                    
                    DAG(A,B)=0;  DAG(B,A)=0; 
                    pdag(A,B)=0;  pdag(B,A)=0;
                    G(A,B)=0;     G(B,A)=0;
                    
                    all_PC{A}=mysetdiff(all_PC{A},B);                    
                    all_can_spouse{A}=myunion(all_can_spouse{A},B);
                    
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
    
    
    % Find V-structures
    for i=1:length(all_can_spouse{A})
        for j=1:length(all_PC{A})
            
            C=all_can_spouse{A}(i);
            B=all_PC{A}(j);
            
            % A->B<-C
                       
            if ~ismember(B,all_sepset{A}{C})
                
                DAG(A,B)=1;     DAG(B,A)=1;
                pdag(A,B) = -1; pdag(C,B) = -1; pdag(B,A) = 0; pdag(B,C) = 0;
                G(A,B) = 1;     G(C,B) = 1;     G(B,A) = 0;    G(B,C) = 0;
                
            end
            
        end
    end
        
    
    [DAG,pdag,G]=meeks(DAG,pdag,G,p);
    
    num_calculated=num_calculated+1;
    
    if num_calculated > length(all_MB{target})
        if ~ismember(1,pdag(target,:))&&~ismember(1,pdag(:,target))
            % P and C have been distinguished
            break;
        end
    end
    
end


P=find(pdag(:,target)==-1);
C=find(pdag(target,:)==-1);
UN=find(pdag(target,:)==1);

PC=myunion(myunion(P,C),UN);

time=toc(start);


