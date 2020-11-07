
function [PC,test,time,sepset] = GetPC_Z(Data,target,alpha,already_calculated_PCD,all_PCD,smaples,p,maxK)
%
% GetPC_Z finds the parents and children of target node on continuous data
%
% INPUT :
%       Data is the data matrix
%       target is the index of target node
%       alpha is the significance level
%       already_calculated_PCD records all nodes that have found PC
%       all_PCD stores the PCDs of all nodes
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


if (nargin == 5)
   [smaples,p]=size(Data);
   maxK=3;
end

start=tic;
test = 0;

PC=[];

if already_calculated_PCD(target)
    [PCD,GetPCD_ntest,~,sepset]=GetPCD_Z(Data,target,alpha, smaples, p, maxK);
    test = test + GetPCD_ntest;
    already_calculated_PCD(target)=0;
    all_PCD{target}=PCD;
else
    PCD=all_PCD{target};
end



for i=1:length(PCD)
    X = PCD(i);
    
    if already_calculated_PCD(X)
        [PCX,GetPCDX_ntest]=GetPCD_Z(Data,X,alpha, smaples, p, maxK);
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





