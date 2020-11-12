
function [MB,P,C,PC,SP,p,SP_sub] = STA(DAG)
%
% STA finds the Markov blankets, parents, children, the union of parents 
% and children, and spouses of each node according to the true network
% 
% INPUT :
%       DAG is the true network (an adjacency matrix)
%
% OUTPUT:
%       MB is the Markov blanket of each node
%       P is the parents of each node
%       C is the children of each node
%       PC is the union of the parents and children of each node
%       SP is the spouses of each node
%       p is the number of nodes
%       SP_sub is the spouses of each node with regard to its common child
%
%

[~,p]=size(DAG);
MB = cell(1,p);
P = cell(1,p);
C = cell(1,p);
PC = cell(1,p);
SP = cell(1,p);

SP_sub=cell(p,p);

for i=1:p
    
    sp_tmp=[];
    p_tmp = parents(DAG, i);
    c_tmp = children(DAG, i, 1);
    pc_tmp = myunion(p_tmp,c_tmp);
    
    P{i}=p_tmp;
    C{i}=c_tmp;
    PC{i}=pc_tmp;
    
    for j=1:length(c_tmp)
        X = c_tmp(j);
        sp_tmp = myunion(sp_tmp,parents(DAG, X));
        
        SP_sub{i,X}=mysetdiff(parents(DAG, X),i);
    end
    sp_tmp = mysetdiff(sp_tmp,i);
    SP{i}=sp_tmp;
    
    mb_tmp= myunion(pc_tmp,sp_tmp);
    MB{i}=mb_tmp;
    
%     i
%     p_tmp
%     c_tmp
%     pc_tmp
%     sp_tmp
%     mb_tmp
%     
end


