function [arrhd_F1,arrhd_precision,arrhd_recall,SHD,reverse,miss,extra]=evaluation_GCS(dag,G)



[SHD,reverse,miss,extra,undirected]=eva_GCS_SHD(G,dag);

[arrhd_F1,arrhd_precision,arrhd_recall]=eva_GCS_arrhd(G,dag,undirected,reverse,miss,extra,SHD);

