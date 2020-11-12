
function [IDT,idT3,idT3_count,test,already_calculated_MB,pc,time] = CMB_subroutine_G2(data,ns,target,alpha,IDT,maxK,already_calculated_MB,all_MB)


start=tic;

test=0;

num=1000;
idT3=zeros(num,2);
idT3_count=0;
idT4=zeros(num,2);
idT4_count=0;

[~,p]=size(data);

[pc,ntest_PC] = HITONPC_G2 (data,target,alpha, ns, p, maxK);

test=test+ntest_PC;

Z=[];

[IDT,idT3,idT4,ntest_CS1,idT3_count,idT4_count]=CausalSearch_G2(data,ns,target,alpha,pc,Z,IDT,idT3,idT3_count,idT4,idT4_count);
test=test+ntest_CS1;

% Step 2: Further test variables with idT = 4

for i=1:idT4_count
    X=idT4(i,1);
    Y=idT4(i,2);

    % optimization£º Do not repeat   2019/3/4
    
    if already_calculated_MB(X)
        
        [all_MB{X},ntest_MB] = HITONMB_G2 (data,X,alpha, ns, p, maxK);
        
        test=test+ntest_MB;
        
        already_calculated_MB(X)=0;
    end
    
    Z=mysetdiff(mysetdiff(all_MB{X},target),Y);
    
    [m,n]=size(Z);
    if m>1&&n==1
        Z=Z';           % Sometimes, the dimension of matrix Z is m*1
    end

    [IDT,idT3,idT4,ntest_CS2]=CausalSearch_G2(data,ns,target,alpha,pc,Z,IDT,idT3,idT3_count,idT4,idT4_count);
    test=test+ntest_CS2;
    
    % no element of IDT is equal to 4
    if isempty(find(IDT==4, 1))
        break;
    end
end

% LINE 12-14 of pseudo code

% for i=1:idT4_count
%     if IDT(idT4(idT4_count,1))==1
%         new_X=idT4(idT4_count,1);
%         new_Z=idT4(idT4_count,2);
%         
%         new_Y_counts=find(idT4(:,2)==new_Z);
%         for j=1:length(new_Y_counts)
%             new_Y_count=new_Y_counts(j);
%             new_Y=idT4(new_Y_count,1);
%             if IDT(Y)==1 && new_Y~=new_X
%                 IDT(new_Z)=1;
%             end
%         end
%              
%     end
% end


parents=find(IDT==1);

for i=1:(length(parents)-1)
    X=parents(i);
    for j=(i+1):length(parents)

        Y=parents(j);
        for k=1:idT4_count
            if idT4(k,1)==X
                Z=idT4(k,2);
                                
                for l=1:idT4_count
                    
                    if k==l
                        continue;
                    end
                    
                    if (idT4(l,2)==Z&&idT4(l,1)==Y)||(idT4(l,1)==Z&&idT4(l,2)==Y)
                        IDT(Z)=1;
                    end
                end
            elseif idT4(k,2)==X
                Z=idT4(k,1);
                
                for l=1:idT4_count
                    
                    if k==l
                        continue;
                    end
                    
                    if (idT4(l,2)==Z&&idT4(l,1)==Y)||(idT4(l,1)==Z&&idT4(l,2)==Y)
                        IDT(Z)=1;
                    end
                end
            end
        end
    end
end


% LINE 15 of pseudo code

IDT(IDT==4)=3;


time=toc(start);





