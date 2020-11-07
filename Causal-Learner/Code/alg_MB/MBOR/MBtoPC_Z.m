function [PC,test,time] = MBtoPC_Z(Data,target,alpha,samples,p,maxK)

start=tic;

test=0;

[MB,ntest_CorrectMB]=CorrectMB_Z(Data,target,alpha, samples, p, maxK);

test = test+ntest_CorrectMB;
PC=MB;

for i=length(MB):-1:1
    X = MB(i);
    MBX = mysetdiff(MB,X);
    cutSetSize = 0;
    
    break_flag=0;
    while length(MBX) >= cutSetSize  && cutSetSize<=maxK
        SS = subsets1(MBX, cutSetSize);  
        for si=1:length(SS)
            Z = SS{si};
            test=test+1;
            [CI]=my_fisherz_test(X,target,Z,Data,samples,alpha);      
            if CI
                PC = mysetdiff(PC,X);
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

PC = sort(PC);

time=toc(start);





