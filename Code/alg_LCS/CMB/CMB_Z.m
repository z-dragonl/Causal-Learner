
function [P,C,PC,UN,test,time] = CMB_Z(Data,target,alpha,samples,p,maxK)
%
% CMB_Z finds and distinguishes the parents and children of target node on continuous data
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

Tmp=[];
Q=target;

DAG=zeros(p,p);
pdag=zeros(p,p);
G=zeros(p,p);

all_idT3=cell(1,p);
all_idT3_count=zeros(1,p);

already_calculated=ones(1,p);

already_calculated_MB=ones(1,p);
all_MB=cell(1,p);


% Step 1: Establish initial ID 
IDT = zeros(p, p);

% if no element of IDT is equal to 3, break;
while length(Tmp)<=p && ~isempty(Q)


    A=Q(1);    Q(1)=[];
    
    if ismember(A,Tmp)
        continue;
    else
        Tmp=[Tmp A];
    end
    
    if already_calculated(A)
        [IDT(A,:),all_idT3{A},all_idT3_count(A),ntest_CMBs1,already_calculated_MB] = CMB_subroutine_Z(Data,samples,A,alpha,IDT(A,:),maxK,already_calculated_MB,all_MB);
        already_calculated(A)=0;
        
        test=test+ntest_CMBs1;
    end

    IDT_A_3=find(IDT(A,:)==3);
    IDT_A_2=find(IDT(A,:)==2);
    IDT_A_1=find(IDT(A,:)==1);
    
    DAG(A,IDT_A_3)=1;  DAG(IDT_A_3,A)=1;
    DAG(A,IDT_A_2)=1;  DAG(IDT_A_2,A)=1;
    DAG(A,IDT_A_1)=1;  DAG(IDT_A_1,A)=1;
    
    pdag(A,IDT_A_3)=1;   pdag(IDT_A_3,A)=1;
    pdag(A,IDT_A_2)=-1;  pdag(IDT_A_2,A)=0;
    pdag(A,IDT_A_1)=0;   pdag(IDT_A_1,A)=-1;
    
    G(A,IDT_A_3)=1;  G(IDT_A_3,A)=1;
    G(A,IDT_A_2)=1;  G(IDT_A_2,A)=0;
    G(A,IDT_A_1)=0;  G(IDT_A_1,A)=1;
    
    

    if ~ismember(1,pdag(target,:))&&~ismember(1,pdag(:,target))
        % P and C have been distinguished
        break;
    end

    % Step 3: Resolve variable set with idT = 3
    
    IDT3_count=find(IDT(A,:)==3);

    for i=1:length(IDT3_count)
        X=IDT3_count(i);
        
        Q=[Q X];  
        
        if already_calculated(X)
            [IDT(X,:),all_idT3{X},all_idT3_count(X),ntest_CMBs2,already_calculated_MB] = CMB_subroutine_Z(Data,samples,X,alpha,IDT(X,:),maxK,already_calculated_MB,all_MB);
            already_calculated(X)=0;
            
            test=test+ntest_CMBs2;
        end
        
        % update IDT according to IDX
        if IDT(X,A)==2
            IDT(A,X)=1;
            for j=1:all_idT3_count(X)
                if all_idT3{X}(j,1)==X
                    Y=all_idT3{X}(j,2);
                    IDT(A,Y)=2;
                elseif all_idT3{X}(j,2)==X
                    Y=all_idT3{X}(j,1);
                    IDT(A,Y)=2;
                end
            end
        end
        
        IDT_X_3=find(IDT(X,:)==3);
        IDT_X_2=find(IDT(X,:)==2);
        IDT_X_1=find(IDT(X,:)==1);

        DAG(X,IDT_X_3)=1;  DAG(IDT_X_3,X)=1;
        DAG(X,IDT_X_2)=1;  DAG(IDT_X_2,X)=1;
        DAG(X,IDT_X_1)=1;  DAG(IDT_X_1,X)=1;
        
        pdag(X,IDT_X_3)=1;   pdag(IDT_X_3,X)=1;
        pdag(X,IDT_X_2)=-1;  pdag(IDT_X_2,X)=0;
        pdag(X,IDT_X_1)=0;   pdag(IDT_X_1,X)=-1;
        
        G(X,IDT_X_3)=1;  G(IDT_X_3,X)=1;
        G(X,IDT_X_2)=1;  G(IDT_X_2,X)=0;
        G(X,IDT_X_1)=0;  G(IDT_X_1,X)=1;
        

              
    end
    
    
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




