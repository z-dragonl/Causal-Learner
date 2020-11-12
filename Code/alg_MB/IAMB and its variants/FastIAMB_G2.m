
function [MB,test,time] = FastIAMB_G2(Data,target,alpha,ns,p,maxK)
%
% FastIAMB_G2 finds the Markov blanket of target node on discrete data
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

MB=[];

test=0;

mb_length_tmp = -1;
tmp_CanMB=[];
CanMB = mysetdiff([1:p],target);
S=[];
pval=zeros(1,p);
dep=zeros(1,p);

for i=1:length(CanMB)
    X = CanMB(i);
    [pval(X),dep(X)]=my_g2_test(X,target,[],Data,ns,alpha);
    test=test+1;
    if isnan(pval(X)) || pval(X)<=alpha
        S = [S,X];
    end
end

TMP=-1;

while( ~isempty(S) )
    
    insufficient_data=0;
    [dep_sort,var_index]=sort(dep(S),'descend');   

    %---------------------------------------------------------
    %  Growing phase.
    
    for i=1:length(S)
        X = S(var_index(i));
        N = size(Data,1);
        qi=ns(MB);
        tmp=[1 cumprod(qi(1:end-1))];
        qs=1+(qi-1)*tmp';
        if isempty(qs),
            df=prod(ns([X target])-1)*prod(ns(MB));
        else
            %   Addition ends
            df=prod(ns([X target])-1)*qs;
        end
        
        if (N>=5*df)
            MB = [MB,X];
        else
            % Not enough data to perform the test
            insufficient_data=1;
            break;
        end
    end
    
    % -----------------------------------------------------------
    % Shrinking phase
    
    remove_flag = 1;
    tmp_MB = MB;
    for i=1:length(MB)
        condset=mysetdiff(MB,MB(i));
        [pval]=my_g2_test(MB(i),target,condset,Data,ns,alpha);
        test=test+1;
        if pval>alpha
            tmp_MB=mysetdiff(tmp_MB,MB(i));
            remove_flag=0;
        end
    end
    
    MB = tmp_MB;
    
    if isequal(TMP,MB)
        break;         
    end
    TMP=MB;
    
    % -----------------------------------------------------------
    % Get new S set
    
    if insufficient_data && remove_flag
        break;
    else
        pval=zeros(1,p);
        dep=zeros(1,p);
        S=[];
        UTBT = mysetdiff(CanMB,MB);
        for i=1:length(UTBT)
            X = UTBT(i);
            [pval(X),dep(X)]=my_g2_test(X,target,MB,Data,ns,alpha);
            test=test+1;
            if isnan(pval(X)) || pval(X)<=alpha
                S = [S,X];
            end
        end
    end

end

time=toc(start);





