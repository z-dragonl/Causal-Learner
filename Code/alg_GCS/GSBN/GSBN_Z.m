

function [G,time,test] = GSBN_Z(Data,alpha,samples,p,maxK)
%
% GSBN_Z learns a causal graph on continuous data
%
% INPUT :
%       Data is the data matrix
%       alpha is the significance level
%       samples is the number of data samples
%       p is the number of nodes
%       maxK is the maximum size of conditioning set
%
% OUTPUT:
%       G is the causal graph
%       time is the runtime of the algorithm
%       test is the number of conditional independence tests
%
%


if (nargin == 2)
   [samples,p]=size(Data);
   maxK=3;
end

start=tic;
test=0;

all_mb=cell(1,p);
all_Bd=cell(1,p);

PP=zeros(p,p);


DAG=zeros(p,p);
pdag=zeros(p,p);
G=zeros(p,p);

for i=1:p
     
    MB=GS_Z(Data,i,alpha,samples,p,maxK);
    PP(i,MB)=1;

end


% AND rule for the correctness of MB

PP_AND = PP;
PP_AND((PP~=PP'))=0;



for i=1:p
    all_mb{i}=find(PP_AND(i,:)==1);
end


% all_mb

% removes the possible spouse links between linked variables X and Y

for i=1:p
    X=i;
    for j=1:length(all_mb{X})
        Y=all_mb{X}(j);
        
        B=myunion(mysetdiff(all_mb{X},Y),mysetdiff(all_mb{Y},X));   
        
        break_flag=0;
        cutSetSize = 0;
        while length(B) >= cutSetSize&&cutSetSize<=maxK
            SS = subsets1(B, cutSetSize);  
            for si=1:length(SS)
                S = SS{si};
                test=test+1;
                [CI]=my_fisherz_test(X,Y,S,Data,samples,alpha);      
                if isnan(CI)       
                    CI=0;
                end
                if(CI==1)            
                    PP_AND(X,Y)=0;  PP_AND(Y,X)=0;
                    break_flag=1;
                    break;
                end
            end
            if( break_flag==1 )
                break;
            end
            cutSetSize = cutSetSize + 1;
        end
        
    end
end


for i=1:p
    all_Bd{i}=find(PP_AND(i,:)==1);
end




DAG=PP_AND;
pdag=PP_AND;
G=PP_AND;


% orient egdes

for i=1:p
    X=i;
    for j=1:length(all_Bd{X})
        Y=all_Bd{X}(j);
        
        SZ=mysetdiff(mysetdiff(all_Bd{X},all_Bd{Y}),Y);
        for k=1:length(SZ)
            Z=SZ(k);
            
            % orient Y -> X
            PP_AND(Y,X)=-1;%PP_AND(X,Y)=0;
            
           
            
            B=myunion(mysetdiff(all_mb{Y},Z),mysetdiff(all_mb{Z},Y));  
            
            break_flag=0;
            cutSetSize = 0;
            while length(B) >= cutSetSize&&cutSetSize<=maxK
                SS = subsets1(B, cutSetSize);   
                for si=1:length(SS)
                    S = SS{si};
                    test=test+1;
                    [CI]=my_fisherz_test(Y,Z,myunion(S,X),Data,samples,alpha);       
                    if isnan(CI)       
                        CI=0;
                    end
                    if(CI==1)            
                        PP_AND(Y,X)=1; %PP_AND(X,Y)=1;
                        break_flag=1;
                        break;
                    end
                end
                if( break_flag==1 )
                    break;
                end
                cutSetSize = cutSetSize + 1;
            end
            
            if  PP_AND(Y,X)==-1
                pdag(Y,X)=-1;pdag(X,Y)=0;    G(Y,X)=1;G(X,Y)=0;
                break;
            end
                        
        end
        
    end
end

[DAG,pdag,G]=meeks(DAG,pdag,G,p);

% draw_graph(G);

G=cpdag_to_dag(G);

time=toc(start);
