function [PC,test,time] = MBtoPC_G2(Data,target,alpha,ns,p,maxK)

start=tic;

test=0;

[MB,ntest_CorrectMB]=CorrectMB_G2(Data,target,alpha, ns, p, maxK);

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
            [pval]=my_g2_test(X,target,Z,Data,ns,alpha);        
            if pval>alpha
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





