
function [PC,test,time,sepset] = GetPC_G2(Data,target,alpha,already_calculated_PCD,all_PCD,ns,p,maxK)
%
% GetPC_G2 finds the parents and children of target node on discrete data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       already_calculated_PCD records all nodes that have found PC
%       all_PCD stores the PCDs of all nodes
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


if (nargin == 5)
   ns=max(Data);
   [~,p]=size(Data);
   maxK=3;
end

start=tic;
test = 0;

PC=[];

if already_calculated_PCD(target)
    [PCD,GetPCD_ntest,~,sepset]=GetPCD_G2(Data,target,alpha, ns, p, maxK);
    test = test + GetPCD_ntest;
    already_calculated_PCD(target)=0;
    all_PCD{target}=PCD;
else
    PCD=all_PCD{target};
end



for i=1:length(PCD)
    X = PCD(i);
    
    if already_calculated_PCD(X)
        [PCX,GetPCDX_ntest]=GetPCD_G2(Data,X,alpha, ns, p, maxK);
        test = test + GetPCDX_ntest;
        already_calculated_PCD(X)=0;
        all_PCD{X}=PCX;
    else
        PCX=all_PCD{X};
    end
    
    if find( PCX== target )
        PC = myunion(PC,X);
    end
    
end

time=toc(start);





