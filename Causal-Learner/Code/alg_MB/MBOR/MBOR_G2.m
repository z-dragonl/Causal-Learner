
function [MB,test,time] = MBOR_G2(Data,target,alpha,ns,p,maxK)
%
% MBOR_G2 finds the Markov blanket of target node on discrete data
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

% ------------------------------------------------------------
% Phase I: Find MB superset (MBS)

[PCS,ntest_PCS,dSep]=PCSuperSet_G2(Data, target, alpha, ns, p, maxK);

[SPS,ntest_SPS] = SPSuperSet_G2(Data, target, alpha, PCS, dSep, ns, p, maxK);

test=test+ntest_PCS;
test=test+ntest_SPS;
MBS=myunion(PCS,SPS);

% D = D(MBS ¡È T)

MBST = myunion(MBS,target);
data=Data(:,MBST);


% ------------------------------------------------------------
% Phase II: Find parents and children of the target


new_ns=ns(:,MBST);
new_p=length(MBST);

target_index = find( MBST==target );
[PC_index,ntest_PC]=MBtoPC_G2(data,target_index, alpha, new_ns, new_p, maxK);       % use new data

PC=MBST(PC_index);

test=test+ntest_PC;
PCSPC = mysetdiff(PCS,PC);

for i=1:length(PCSPC)
    X = PCSPC(i);
    
    X_index = find( MBST==X );
    [PCX_index,ntest_PCX]=MBtoPC_G2(data,X_index, alpha, new_ns, new_p, maxK);      % use new data
    
    PCX=MBST(PCX_index);
    
    test=test+ntest_PCX;
    if ismember(target,PCX)
        PC = myunion(PC,X);
    end
end

% ------------------------------------------------------------
% Phase III: Find spouses of the target

SP = [];
for i=1:length(PC)
    X = PC(i);
    [PCX_SP,ntest_PCX_SP]=MBtoPC_G2(Data,X, alpha, ns, p, maxK);
    test=test+ntest_PCX_SP;
    XDPCT = mysetdiff( PCX_SP,myunion(PC,target) );
    for j=1:length(XDPCT)
        Y = XDPCT(j);
        MBSTY = mysetdiff( MBS,myunion(target,Y) );
        
        tmp_dep = 10000;
        tmp_Z = [];
        cutSetSize = 0;
        while length(MBSTY) >= cutSetSize && cutSetSize<=maxK
            SS = subsets1(MBSTY, cutSetSize);    
            for si=1:length(SS)
                Z = SS{si};                
                test=test+1;
                [pval,dep]=my_g2_test(Y,target,Z,Data,ns,alpha);        
                if pval>alpha
                    if dep<tmp_dep
                        tmp_dep = dep;
                        tmp_Z = Z;
                    end
                end
            end
            cutSetSize = cutSetSize + 1;
        end
        if tmp_dep~=10000
            test=test+1;
            [pval_SP]=my_g2_test(Y,target,myunion(tmp_Z,X),Data,ns,alpha);
            if pval_SP<=alpha||isnan(pval_SP)
                SP = myunion(SP,Y);
            end
        end
        
    end
end

MB=myunion(PC,SP);

time=toc(start);









