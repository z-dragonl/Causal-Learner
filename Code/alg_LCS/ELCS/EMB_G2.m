function [MB, pc, parent, child, unsure, SpouseY, ntest, ntime] = EMB_G2(Data, target, alpha, ns, p, maxk)

% Input :
%       Data is the data matrix, NbVar columns * N rows
%       target is the index of target variable
%       alpha is the significance level
%       ns = max(data)
%       p is the number of variables
%       maxk is the caidinal of max subset
%output:
%       MB is the parents, children and spouses of the target
%       pc is the parents and children of the target
%       parent is the parents of the target
%       children is the children of the target
%       unsure is the undistinguished PC variables of the target
%       SpouseY is the set that contains the  spouse sets of all PC variables 
%       ntest is the number of conditional independence
%       time is the runtime of the algorithm

    start = tic;
    ntest = 0;
    Spouse = [];

    %------------------Step1: use other algorithms to find pc of variable target-------------------------% 
    [pc, ntesthpc, ~,sepset] = HITONPC_G2(Data,target,alpha, ns, p, maxk);
    ntest = ntest + ntesthpc;

    %------------------Step2-I: find the candidate spouse of each pci in pc -------------------------------% 
    SpouseY = cell(1,p);
    CanSp = mysetdiff(1:p, [pc,target]);
%     mysepsetpc = cell(1,p);
    for i = 1 : length(CanSp)
        tmp = [];
        for j = 1 : length(pc)
            ntest = ntest + 1;
            [pvaltmp] = my_g2_test(CanSp(i), pc(j), [], Data, ns,alpha);        % use G^2 to test conditional independence
            if isnan(pvaltmp)
                CI = 0;
            else
               if pvaltmp <= alpha
                  CI = 0;
               else
                  CI = 1;
               end
            end
            if CI == 0
                tmp = [tmp, pc(j)];
            end
        end

        ntest = ntest + 1;
        [pvaltmp] = my_g2_test(CanSp(i), target, tmp, Data, ns,alpha);        % use G^2 to test conditional independence
        if isnan(pvaltmp)
            CI = 0;
        else
           if pvaltmp <= alpha
              CI = 0;
           else
              CI = 1;
           end
        end
        if CI == 1
            continue;
        else
            rest  = tmp;
            for k = 1 : length(rest)
                ntest=ntest+1;
                [pval]=my_g2_test(target,CanSp(i),[sepset{1,CanSp(i)},rest(k)],Data,ns,alpha);        % use G^2 to test conditional independence
                if isnan(pval)
                    CI=0;
                else
                   if pval<=alpha
                      CI=0;
                   else
                      CI=1;
                   end
                end
                if CI==0
                    SpouseY{1,rest(k)} = myunion(SpouseY{1,rest(k)}, CanSp(i));
                end
            end
        end
    end

    SY = cell(1,length(pc));
    for i = 1 : length(pc)
        SY{1,i} = SpouseY{1,pc(i)};
    end
    SpouseY = SY;

    SpouseYCheckC = SY; %used to check some children 
    %------------------Step2-II: remove false spouse of each pci in pc -------------------------------------% 
    %remove  flase spouse of children (intra )
    for i = 1 : length(pc)
        SpYtmp = SpouseY{1, i};
        %-----------------sort: to reduce the time------------------
        DepSpY = zeros(1,p);
        for j = 1 : length(SpYtmp)
            ntest = ntest + 1;
            [pval,DepSpY(SpYtmp(j))]=my_g2_test(SpYtmp(j), pc(i), [], Data, ns,alpha);        % use G^2 to test conditional independence
        end
    %     
        [dep_sort,var_index]=sort(DepSpY(SpYtmp),'descend');

        CanSptmp = [];
        for j = 1 : length(SpYtmp)
            CanSptmp = [CanSptmp, SpYtmp(var_index(j))];
        end
        SpYtmp = CanSptmp;
        SpY = [];

        %-----------------------------------------------------------
        for m=1:length(SpYtmp)
            Y=SpYtmp(m);
            SpY=[SpY Y];
            pc_length=length(SpY);
            pc_tmp=SpY;
            for j=pc_length:-1:1
                X=SpY(j);
                CanPC=mysetdiff(pc_tmp, X);
                break_flag=0;
                cutSetSize = 0;
                while length(CanPC) >= cutSetSize &&cutSetSize<=maxk

                    SS = subsets1(CanPC, cutSetSize);    % all subsets of size cutSetSize
                    for si=1:length(SS)
                        Z = SS{si};

                        if X~=Y             
                            if isempty(find(Z==Y, 1))
                                continue;
                            end
                        end  

                        ntest=ntest+1;
                        [pval]=my_g2_test(X,pc(i),Z,Data,ns,alpha);        % use G^2 to test conditional independence
                        if isnan(pval)
                            CI=0;
                        else
                           if pval<=alpha
                              CI=0;
                           else
                              CI=1;
                           end
                        end
                        if CI==1
                            pc_tmp=CanPC;
        %                     sepset{1,X}=Z;
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
            SpY=pc_tmp;
        end
        SpouseY{1, i} = SpY;
    end 

    for i = 1 : length(pc)
        X = pc(i);
        SpYtmp = SpouseY{1, i};
        for j = 1 : length(SpYtmp)
            Y = SpYtmp(j);
            pcRest = mysetdiff(pc, pc(i));
            ConSet = [pcRest, target, mysetdiff(SpYtmp, SpYtmp(j))];
            CanPC=ConSet;
            break_flag=0;
            cutSetSize = 0;
            while length(CanPC) >= cutSetSize &&cutSetSize<=maxk

                SS = subsets1(CanPC, cutSetSize);    % all subsets of size cutSetSize
                for si=1:length(SS)
                    Z = SS{si};
                    ntest=ntest+1;
                    [pval]=my_g2_test(pc(i), SpYtmp(j),Z,Data,ns,alpha);        % use G^2 to test conditional independence
                    if isnan(pval)
                        CI=0;
                    else
                       if pval<=alpha
                          CI=0;
                       else
                          CI=1;
                       end
                    end
                    if CI==1
                        SpouseY{1,i} = mysetdiff(SpouseY{1,i},SpYtmp(j));
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

    %------------------Step2-III: remove false spouse of parent (  condition  Ind error)---------------------% 
    % remove  false spouseS of parent 
    for i = 1 : length(pc)
        SpYtmp  = SpouseY{1, i};
        pcRest = mysetdiff(pc, pc(i));
        CanSpY = SpouseY{1, i}; 
        for j = 1 : length(SpYtmp)
            cutSetSize = 0;    
            while length(pcRest)>=cutSetSize &&cutSetSize<=maxk
                flags = 1;
                SS = subsets1(pcRest, cutSetSize);    % all subsets of size cutSetSize
                for si = 1 : length(SS)
                    ConSet = [SS{si}, pc(i)];
                    ntest = ntest + 1;
                    [pval]=my_g2_test(target, SpYtmp(j), ConSet, Data, ns,alpha);        % use G^2 to test conditional independence
                    if isnan(pval)
                        CI = 0;
                    else
                       if pval <= alpha
                          CI = 0;
                       else
                          CI = 1;
                          flags = 0;
                          CanSpY = mysetdiff(CanSpY, SpYtmp(j));
                          break;
                       end
                    end
                end
                if flags == 0
                    break;
                end
                cutSetSize = cutSetSize + 1;
            end
        end
        SpouseY{1, i} = CanSpY;
    end


    %------------------Step3: remove false pc---------------------------- -------------------------------% 
    CanPC = pc;
    for i = 1 : length(pc)
        CanSpY = SpouseY{1, i};
        pcRest = mysetdiff(pc, pc(i));
        cutSetSize = 1;    
        while length(pcRest)>=cutSetSize &&cutSetSize<=maxk
            flags = 1;
            SS = subsets1(pcRest, cutSetSize);    % all subsets of size cutSetSize
            for si = 1 : length(SS)
                ConSet = [SS{si}, CanSpY];
                ntest = ntest + 1;
                [pval]=my_g2_test(target, pc(i), ConSet, Data, ns,alpha);        % use G^2 to test conditional independence
                if isnan(pval)
                    CI = 0;
                else
                   if pval <= alpha
                      CI = 0;
                   else
                      CI = 1;
                      flags = 0;
                      CanPC = mysetdiff(CanPC, pc(i));
                      SpouseY{1, i} = [];
                      break;
                   end
                end
            end
            if flags == 0
                break;
            end
            cutSetSize = cutSetSize + 1;
        end
    end

    %------------------Step4: DistinguishPC---------------------------- -------------------------------% 
    Spouse = [];
    for i = 1 : length(pc)
        CanSpY = SpouseY{1, i};
        Spouse = myunion(CanSpY, Spouse);
    end

    parent = [];
    child = [];
    unsure = [];
    for i = 1 : length(pc)
        if length(SpouseY{1, i})>0
            child = [child, pc(i)];
        end
    end

    %--------------------------------------------------------------------------------------
    % using N-structures to identify the direct effects  of a given variable
    child_two =[];%Satisfy certain conditions
    for i = 1 : length(pc)
        if length(SpouseY{1, i})==0
            if(intersect(SpouseYCheckC{1,i},Spouse))
                child_two = [child_two,pc(i)];
            end
        end
    end
    %--------------------------------------------------------------------------------------

    child_two = unique(child_two);
    child = [child, child_two];
    child = unique(child);
    parent = [];

    unsure = mysetdiff(CanPC, child);

    % using Lemma 1 (a) to identify the direct causes  of the target  variable
    for i = 1 : length(unsure)
        for j = i+1 : length(unsure)
            ntest = ntest + 1;
            [pval]=my_g2_test(unsure(i), unsure(j), [], Data, ns,alpha);        
            if isnan(pval)
                CI = 0;
            else
               if pval <= alpha
                    CI = 0;
               else
                    CI = 1;
                    ntest = ntest + 1;
                    [pval]=my_g2_test(unsure(i), unsure(j), target, Data, ns,alpha);       
                    if isnan(pval)
                        CII = 0;
                    else
                       if pval <= alpha
                          CII = 0;
                       else
                          CII = 1;
                       end
                    end
                    if CII == 0 %conditional dependence 
                        parent = [parent,[unsure(i), unsure(j)]];
                    end
               end
            end
        end
    end
    parent = unique(parent);

    % using Lemma 1 (b) to identify the direct effects  of the target  variable
    rest = mysetdiff(CanPC, [child,parent]);
    if ~isempty(parent)
        for i = 1 : length(rest)
            ntest = ntest + 1;
            [pval]=my_g2_test(parent(1), rest(i), [], Data, ns,alpha);        
            if isnan(pval)
                CI = 0;
            else
               if pval <= alpha
                    CI = 0;
               else
                    CI = 1;
               end
            end
            if CI==0
                ntest = ntest + 1;
                [pval]=my_g2_test(parent(1), rest(i), target, Data, ns,alpha);       
                if isnan(pval)
                    CII = 0;
                else
                   if pval <= alpha
                      CII = 0;
                   else
                      CII = 1;
                   end
                end
                if CII == 1 %conditional dependence 
                    child = [child,rest(i)];
                end
            end
        end
    end
    unsure = mysetdiff(CanPC, [child,parent]);

    pc = CanPC;
    MB = [pc, Spouse];
    MB = sort(MB);
    ntime = toc(start);
end
