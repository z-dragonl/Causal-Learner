function [PCS,test,dSep,time] = PCSuperSet_G2(Data,target,alpha,ns,p,maxK)

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
    [pval]=my_g2_test(X,target,[],Data,ns,alpha);
    if pval>alpha
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
        [pval]=my_g2_test(X,target,Y,Data,ns,alpha);
        if pval>alpha
            tmp_PCS = mysetdiff(tmp_PCS,X);
            dSep{X} = Y;
            
            break;     
        end
    end
end
PCS = tmp_PCS;
time=toc(start);





