function [arrhd_F1,arrhd_precision,arrhd_recall,SHD,reverse,miss,extra,undirected]=evaluation_LCS(P,C,UN,target,G)

[~,trueP,trueC,truePC] = STA(G);

PC=myunion(myunion(P,C),UN);

[arrhd_F1,arrhd_precision,arrhd_recall]=eva_LCS_arrhd(PC,P,C,trueP{target},trueC{target},truePC{target});

[~,p]=size(G);

local_G=zeros(p,p); DAG=zeros(p,p);

local_G(target,trueC{target})=1;   local_G(trueP{target},target)=1;
DAG(target,C)=1; DAG(P,target)=1;
DAG(target,UN)=1;DAG(UN,target)=1;

[SHD,reverse,miss,extra,undirected]=eva_LCS_SHD(local_G,DAG);

