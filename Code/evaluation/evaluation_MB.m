

function [adj_F1,adj_precision,adj_recall]=evaluation_MB(learned_MB,target,graph)


[trueMB]= STA(graph);

[adj_F1,adj_precision,adj_recall]=eva_MB_adj(learned_MB,trueMB{target});

