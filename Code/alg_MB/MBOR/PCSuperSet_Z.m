function [PCS,test,dSep,time] = PCSuperSet_Z(Data,target,alpha,samples,p,maxK)

start=tic;
test = 0;

dSep=cell(1,p);

% ------------------------------------------------------------
% Phase I: Remove X if T กอ X

PCS = mysetdiff([1:p],target);
tmp_PCS = PCS;
for i=1:length(PCS)
    X = PCS(i);
    test=test+1;
    [CI]=my_fisherz_test(X,target,[],Data,samples,alpha);
    if CI
        tmp_PCS = mysetdiff(tmp_PCS,X);
        dSep{X} = {};
    end
end
PCS = tmp_PCS;

% ------------------------------------------------------------
% Phase II:Remove X if T กอ X|Y

tmp_PCS = PCS;
for i=1:length(PCS)
    X = PCS(i);
    PCSX = mysetdiff( PCS,X );
    for j=1:length(PCSX)
        Y = PCSX(j);
        test=test+1;
        [CI]=my_fisherz_test(X,target,Y,Data,samples,alpha);
        if CI
            tmp_PCS = mysetdiff(tmp_PCS,X);
            dSep{X} = Y;
            
            break;      
        end
    end
end
PCS = tmp_PCS;
time=toc(start);





