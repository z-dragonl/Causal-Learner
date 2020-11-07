function [SPS,test,time] = SPSuperSet_Z(Data,target,alpha,PCS,dSep,samples,p,maxK)

start=tic;
test = 0;


SPS=[];

for i=1:length(PCS)
    X = PCS(i);
    SPSX=[];
    UTPCS = mysetdiff( [1:p],myunion(PCS,target) );
    for j=1:length(UTPCS)
        Y = UTPCS(j);
        test=test+1;
        [CI]=my_fisherz_test(Y,target,myunion(dSep{Y},X),Data,samples,alpha);
        if CI==0 || isnan(CI)
            SPSX = [SPSX,Y];
        end
    end
    
    tmp_SPSX = SPSX;
    for k=1:length(SPSX)
        Y = SPSX(k);
        SPSXY = mysetdiff(SPSX,Y);
        for l = 1:length(SPSXY)
            Z = SPSXY(l);
            test=test+1;
            [CI]=my_fisherz_test(Y,target,myunion(Z,X),Data,samples,alpha);
            if CI
                tmp_SPSX = mysetdiff(tmp_SPSX,Y);
                break;      
            end
        end
    end
    SPSX=tmp_SPSX;   
        
    SPS = myunion(SPS,SPSX);
end

time=toc(start);





