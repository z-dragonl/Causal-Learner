
function [MYparent,Mychild,PC,Undirected,ntest,ntime] = ELCS_G2(Data,target,alpha,ns,p,maxk)
    start = tic;
    ntest = 0;
%     all_sepset = cell(1,p);
    DAG = zeros(p,p);
    pdag = zeros(p,p);
    G = zeros(p,p);
    PC_calculated = ones(1,p);
    Q = target;
    Tmp = [];
    [MB, pc, parent, child, unsure, SpouseY, test, time] = EMB_G2(Data, target, alpha, ns, p, maxk);
    ntest = ntest + test;
    PC = pc;
    MYparent = parent;
    Mychild = child;
    Myunsure = unsure;
    PC_calculated(target) = 0;
    while length(Tmp)<=p && ~isempty(Q)


        A=Q(1);    Q(1)=[];

        if ismember(A,Tmp)
            continue;
        else
            Tmp=[Tmp A];
        end

        if PC_calculated(A)
            [MB, pc, parent, child, unsure, SpouseY, test, time] = EMB_G2(Data, A, alpha, ns, p, maxk);
            PC_calculated(A)=0;
            ntest = ntest + test;
        end

        for j = 1 : length(parent)
             B = parent(j);
             DAG(A,B)=1;  DAG(B,A)=1;

             if pdag(A,B)==0 && pdag(B,A)==0
                pdag(A,B)=1;  pdag(B,A)=1;
                G(A,B)=1;  G(B,A)=1;
            end


            if pdag(B,A) ==1 && pdag(B,A) ~=-1
                pdag(B,A) = -1;  pdag(A,B) = 0; 
                G(B,A) = 1;    G(A,B) = 0;  
            end
        end

        for j = 1 : length(pc)
            B = pc(j);
            DAG(A,B)=1;  DAG(B,A)=1;
            if pdag(A,B)==0 && pdag(B,A)==0
                    pdag(A,B)=1;  pdag(B,A)=1;
                    G(A,B)=1;  G(B,A)=1;
            end

            for k = 1 : length(child)
                 if pdag(A,child(k)) ==1 
                     pdag(A,child(k)) = -1;  pdag(child(k),A) = 0;
                     G(A,child(k)) = 1;        G(child(k),A) = 0;  
                 end
            end
            if length(SpouseY{1,j})>0
                for k = 1 : length(SpouseY{1,j})
                    C = SpouseY{1,j}(k);
                    DAG(C,B)=1;  DAG(B,C)=1;
                    if pdag(C,B)==0 && pdag(B,C)==0
                            pdag(C,B)=1;  pdag(B,C)=1;
                            G(C,B)=1;  G(B,C)=1;
                    end

    %                 if pdag(A,B)==1&&pdag(C,B)==-1||pdag(A,B)==-1&&pdag(C,B)==1||pdag(A,B)==1&&pdag(C,B)==1
                        pdag(A,B) = -1; pdag(C,B) = -1; pdag(B,A) = 0; pdag(B,C) = 0;
                        G(A,B) = 1;     G(C,B) = 1;     G(B,A) = 0;    G(B,C) = 0;
    %                 end
                end
            end
        end

        for j = 1 : length(unsure)
             B=unsure(j);   Q=[Q B];
             DAG(A,B)=1;  DAG(B,A)=1;
            if pdag(A,B)==0 && pdag(B,A)==0
                    pdag(A,B)=1;  pdag(B,A)=1;
                    G(A,B)=1;  G(B,A)=1;
            end

        end
         [DAG,pdag,G]=meeks(DAG,pdag,G,p);
        if ~ismember(1,pdag(target,:))&&~ismember(1,pdag(:,target))
            % P and C have been distinguished
            break;
        end
    end


    Undirected = [];
    for i = 1 :  length(Myunsure)
        if pdag(Myunsure(i),target)==-1
            MYparent = [MYparent,Myunsure(i)];
        elseif pdag(target,Myunsure(i))==-1
            Mychild = [Mychild,Myunsure(i)];
        else
            Undirected = [Undirected,Myunsure(i)];
        end
    end
    ntime=toc(start);
end

