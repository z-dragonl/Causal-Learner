
function [IDT,idT3,idT4,test,idT3_count,idT4_count,time] = CausalSearch_Z(data,samples,target,alpha,PC,Z,IDT,idT3,idT3_count,idT4,idT4_count)


start=tic;

test=0;

% Step 1: Single PC

if length(PC)==1
    IDT(PC)=3;
end

% Step 2: Check C2 & C3


for i=1:(length(PC)-1)
    X=PC(i);
    for j=(i+1):length(PC)
        Y=PC(j);
        
        test=test+2;
        
        CI_XYZ=my_fisherz_test(X,Y,Z,data,samples,alpha);
        CI_XYZ_T=my_fisherz_test(X,Y,myunion(Z,target),data,samples,alpha);
        
        if isnan(CI_XYZ)
            CI_XYZ=0;
        end
        
        if isnan(CI_XYZ_T)
            CI_XYZ_T=0;
        end
        
        if CI_XYZ && CI_XYZ_T==0
                IDT(X)=1;  IDT(Y)=1;
        elseif CI_XYZ==0 && CI_XYZ_T
            if IDT(X)==1
                IDT(Y)=2;
            elseif IDT(Y)~=2
                IDT(Y)=3;
            end
            
            if IDT(Y)==1
                IDT(X)=2;
            elseif IDT(X)~=2
                IDT(X)=3;
            end
            
            % add (X,Y) to pairs with idT=3
            idT3_count=idT3_count+1;
            idT3(idT3_count,1)=X;
            idT3(idT3_count,2)=Y;
        else
            if (IDT(X)==0&&IDT(Y)==0)||(IDT(X)==4&&IDT(Y)==4)  %???
                IDT(X)=4;  IDT(Y)=4;
            end
            
            % add (X,Y) to pairs with idT=4
            idT4_count=idT4_count+1;
            idT4(idT4_count,1)=X;
            idT4(idT4_count,2)=Y;
        end
       
    end
end


% Step 3: identify idT = 3 pairs with known parents

for i=1:length(PC)
    X=PC(i);
    if IDT(X)==1
        
        for j=1:idT3_count
            if idT3(j,1)==X
                Y=idT3(j,2);
                IDT(Y)=2;
            elseif idT3(j,2)==X
                Y=idT3(j,1);
                IDT(Y)=2;
            end
        end
        
    end
end


time=toc(start);





