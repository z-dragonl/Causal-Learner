
function [P,C,PC,UN,test,time] = PCD_by_PCD_Z(Data,target,alpha,samples,p,maxK)
%
% PCD_by_PCD_Z finds and distinguishes the parents and children of target node on continuous data
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

PC_calculated=ones(1,p);

all_PC=cell(1,p);
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
    
    
    % Get PC(A)
    if PC_calculated(A)
        [all_PC{A},ntest1,~,all_sepset{A}]=MMPC_Z(Data,A,alpha,samples,p,maxK);
        test=test+ntest1;
        PC_calculated(A)=0; 
    end

    
    % Get PC of each node in PC(A)
    for i=1:length(all_PC{A})
        
        B=all_PC{A}(i);   Q=[Q B];
        
        if ~ismember(A,all_PC{B})
            continue;
        end
        
        DAG(A,B)=1;  DAG(B,A)=1;        % A is member of PC(B)  (AND rule)
        
        if pdag(A,B)==0 && pdag(B,A)==0
            pdag(A,B)=1;  pdag(B,A)=1;
            G(A,B)=1;     G(B,A)=1;
        end
        
        
        % Get PC of each node in PC(B)
        for j=1:length(all_PC{B})
            
            C=all_PC{B}(j);
            
            if ismember(C,all_PC{A})||C==A    % C may be the PC of A, and C may be A
                continue;
            end

            if DAG(C,B)==1 && DAG(B,C)==1
                if ~ismember(B,all_sepset{A}{C})
                    
                    pdag(A,B) = -1; pdag(C,B) = -1; pdag(B,A) = 0; pdag(B,C) = 0;
                    G(A,B) = 1;     G(C,B) = 1;     G(B,A) = 0;    G(B,C) = 0;
                    
                end
            end
        end
    end
    
    [DAG,pdag,G]=meeks(DAG,pdag,G,p);
    
    num_calculated=num_calculated+1;
    
    if num_calculated > length(all_PC{target})
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


