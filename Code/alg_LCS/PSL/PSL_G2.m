
function [P,C,PC,UN,ntest,time] = PSL_G2(Data,target,alpha,ns,p,k)

depth=1;

start=tic;

ntest=0;


DAG=zeros(p,p);
pdag=zeros(p,p);
G=zeros(p,p);

PC_calculated=ones(1,p);

all_PC=cell(1,p);
all_sepset=cell(1,p);

Q=target;
Tmp=[];


lier_num=1;
lier_nodes=cell(1,p);
lier_nodes{1}=target;

ling=1;



while length(Tmp)<=p && ~isempty(Q)
    
    
    A=Q(1);    Q(1)=[];
    
    if ismember(A,Tmp)
        continue;
    else
        Tmp=[Tmp A];
    end
    

    
    % Get PC(A)
    if PC_calculated(A)
        [all_PC{A},ntest1,ntime1,all_sepset{A}]=PCsimple_G2 (Data,A,alpha,ns,p,k);
        
        ntest=ntest+ntest1;
        PC_calculated(A)=0; 
    end
    

%     all_PC{A}
    
    % Get PC of each node in PC(A)
    for i=1:length(all_PC{A})
        
        B=all_PC{A}(i);   Q=[Q B];
        
        
        % Get PC(B)
        if PC_calculated(B)
            [all_PC{B},ntest2,ntime2,all_sepset{B}]=PCsimple_G2 (Data,B,alpha,ns,p,k);
            
            ntest=ntest+ntest2;
            PC_calculated(B)=0;
        end
        
        
        % symmetry constraint
%         if ~ismember(A,all_PC{B})   % A is member of PC(B)  (AND rule)
%             continue;
%         end

        
        DAG(A,B)=1;  DAG(B,A)=1;       
        
        if pdag(A,B)==0 && pdag(B,A)==0
            pdag(A,B)=1;  pdag(B,A)=1;
            G(A,B)=1;  G(B,A)=1;
        end

        
        % find collider V-structure  B->A<-B_tmp
        
        if i~=length(all_PC{A})
            for long=(i+1):length(all_PC{A})
                B_tmp=all_PC{A}(long);
                

                if B_tmp==B||ismember(B_tmp,all_PC{B})
                    continue;
                end

 
%                 %add 0909
%                 if B_tmp==B||ismember(B_tmp,all_PC{B})||ismember(B,all_PC{B_tmp})
%                     continue;
%                 end
%                 %add0909
                
                % symmetry constraint
%                 if PC_calculated(B_tmp)
%                     [all_PC{B_tmp},ntest3,all_sepset{B_tmp}]=MMPC_optimize (Data,B_tmp,alpha,ns,p,k);
%                     
%                     ntest=ntest+ntest3;
%                     PC_calculated(B_tmp)=0;
%                 end
%                 
%                 
%                 if ~ismember(A,all_PC{B_tmp})   % A is member of PC(B_tmp)  (AND rule)
%                     continue;
%                 end
%                 if ismember(B,all_PC{B_tmp})
%                     continue;
%                 end
                
                
                if ~ismember(A,all_sepset{B}{B_tmp})
                    ntest=ntest+1;
                    pval=my_g2_test(B,B_tmp,myunion(A,all_sepset{B}{B_tmp}),Data,ns,alpha);
                    
                    if isnan(pval)||pval<=alpha
                        
                        
                        if 1%pdag(B,A)==1&&pdag(B_tmp,A)==1   % only direct undirected edges / result is not good
                            pdag(B,A) = -1; pdag(B_tmp,A) = -1; pdag(A,B) = 0; pdag(A,B_tmp) = 0;
                            G(B,A) = 1;     G(B_tmp,A) = 1;     G(A,B) = 0;    G(A,B_tmp) = 0;
                        end
                        
%                         fprintf('\n collider V %.f->%.f-<%.f\n',B,A,B_tmp);
%                         
%                         close;
%                         draw_graph(G);
%                         pause;
                        
                    end
                end
                
            end
        end
        

        
        % find non-collider V-structure  A->B<-C
        
        for j=1:length(all_PC{B})
            
            C=all_PC{B}(j);
            
            if ismember(C,all_PC{A})||C==A
                continue;
            end
            
            if ~ismember(B,all_sepset{A}{C})
                ntest=ntest+1;
                pval=my_g2_test(A,C,myunion(B,all_sepset{A}{C}),Data,ns, alpha);
                
                if isnan(pval)||pval<=alpha
%                     if ismember(C,SP_ling{A,B})
                    
                    if 1%pdag(A,B) ==1&&pdag(C,B) ==1 % only direct undirected edges / result is not good
                        DAG(C,B)=1;  DAG(B,C)=1;
                    
                        pdag(A,B) = -1; pdag(C,B) = -1; pdag(B,A) = 0; pdag(B,C) = 0;
                        G(A,B) = 1;     G(C,B) = 1;     G(B,A) = 0;    G(B,C) = 0;
                        
                    end
                    
%                     fprintf('\n non-collider V %.f->%.f-<%.f\n',A,B,C);
%                     
%                     close;
%                     draw_graph(G);
%                     pause;
                    
                end
            end

        end
    end
    
    [DAG,pdag,G]=meeks(DAG,pdag,G,p);
    
    
    ling=ling-1;
    
%     A
%     close;
%     draw_graph(G);
%     pause;
    

    if ling==0

        up_num=length(lier_nodes{lier_num});
        
        for i=1:up_num
            up_num_PC=lier_nodes{lier_num}(i);           
            lier_nodes{lier_num+1}=myunion(lier_nodes{lier_num+1},all_PC{up_num_PC});
        end
        
        ling=ling+length(mysetdiff(lier_nodes{lier_num+1},Tmp));
        
        lier_num=lier_num+1;

%         close;
%         draw_graph(G);
%         pause;
        
    end


    if lier_num>depth

        break_flag=1;
        for i=1:length(lier_nodes{depth})
            up_A=lier_nodes{depth}(i);
           
            if ismember(1,pdag(up_A,:))
                break_flag=0;
            end
        end
        
        if break_flag
            break;
        end
    end
    

end

G=cpdag_to_dag(G);

time=toc(start);

P=find(pdag(:,target)==-1);
C=find(pdag(target,:)==-1);
UN=find(pdag(target,:)==1);

PC=myunion(myunion(P,C),UN);
% for i=1:(lier_K+2)
%     lier_nodes{i}
% end

close;

