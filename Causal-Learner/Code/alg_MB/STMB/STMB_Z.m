
function [MB,test,time] = STMB_Z(Data,target,alpha,samples,p,maxK)
%
% STMB_Z finds the Markov blanket of target node on continuous data
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
%       MB is the Markov blanket of the target
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

[PC, ntest_PC, time_PC, TSepset] = PCsimple_Z(Data, target, alpha, samples, p, maxK);
% PC

test = test+ntest_PC;
spouse = cell(1,p);
remove = [];

NonTPC = mysetdiff(1:p,myunion(PC,target));

%-------------------------------------------------------- step1
% find spouses and remove non-child descendants

for i = 1:length(PC)
    Y = PC(i);
    break_flag=0;
    for j=1:length(NonTPC)
        X = NonTPC(j);
        
        S =  myunion(TSepset{1,X}, Y);    
        test = test +1;
        [CI]=my_fisherz_test(target,X,S,Data,samples,alpha);     
        if isnan(CI)       
            CI=0;
        end
        
        if(CI==0)     
            SZ = mysetdiff( myunion(PC,X),Y );
            cutSetSize = 0;
            while length(SZ) >= cutSetSize&&cutSetSize<=maxK
                SS = subsets1(SZ, cutSetSize);   
                for si=1:length(SS)
                    S = SS{si};
                    test=test+1;
                    [CI]=my_fisherz_test(Y,target,S,Data,samples,alpha);    
                    if isnan(CI)        
                        CI=0;
                    end
                    if(CI==1)           
                        remove = myunion(remove,Y);
                        spouse{1,Y} = [];
                        break_flag=1; 
                        break;
                    end
                end
                if( break_flag==1 )
                    break;
                end
                cutSetSize = cutSetSize + 1;       
            end
            if( break_flag==1 )
                break;
            end
            
            spouse{1,Y} = myunion ( spouse{1,Y} ,X );       % line 13 of STMB
        end
    end
end
PC = mysetdiff( PC,remove );

%-------------------------------------------------------- step2
% test for false positive spouses X

for a = 1:length(PC)

    i=PC(a);
    spousety = spouse{1,i};
    
    if (~isempty(spousety))
        Y=i;
        for j=1:length(spousety)
            X = spousety(j);
            
            testset = mysetdiff( myunion(PC,spouse{1,Y}),X );
            
            test=test+1;
            [CI]=my_fisherz_test(X,target,testset,Data,samples,alpha);      
            
            if isnan(CI)      
                CI=0;
            end
            
            if(CI==1)            
                spouse{1,Y}=mysetdiff( spouse{1,Y},X );
            end
        end
        
    end
    
end

%-------------------------------------------------------- step3
% test for other non-MB descendants X in the PC set

M = PC;

for i = 1:length(M)
    X = M(i);
    
    pc_tmp=PC;
    CanM = mysetdiff(pc_tmp,X);
    
    spouse_test=[];
    
    for maxK=1:length(CanM)
        pc_var = CanM(maxK);
        spouse_test = myunion(spouse_test,spouse{1,pc_var});
    end
    TestSet = myunion(CanM,spouse_test);
    
    test=test+1;
    [CI]=my_fisherz_test(X,target,TestSet,Data,samples,alpha);     
    if isnan(CI)
        CI=0;
    end
    if CI==1      
        PC=mysetdiff( PC,X );
        spouse{1,X} = [];
    end
end


MB=myunion(PC,cell2mat(spouse));

time=toc(start);


