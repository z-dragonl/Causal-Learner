function [MB,test,time] = EAMB_G2(Data,target,alpha,ns,p,maxK)

k_or=0.05;
test=0;


ns=max(Data);
[~,n_vars]=size(Data);

start=tic;

[MB,ntest,~,index] = SRFS_main_G2(Data,target,alpha, ns, n_vars);
test= test+ntest;

or_len=round(length(index)*k_or);

for i=1:or_len
    [mb,ntest] = SRFS_minor_G2(Data,index(i),alpha, ns, n_vars);
    test = test+ntest;
    if ismember(target,mb)
        MB=[MB index(i)];
    end
end

time=toc(start);

end

