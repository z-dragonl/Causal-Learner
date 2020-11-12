

% ************ Names of the 31 networks are as follows ************ 
   
% ----- Discrete Bayesian Networks

% Small Networks (<20 nodes)   asia, cancer, earthquake, sachs, survey

% Medium Networks (20每50 nodes)   alarm, barley, child, insurance, mildew, water

% Large Networks (50每100 nodes)   hailfinder, hepar2, win95pts

% Very Large Networks (100每1000 nodes)   andes, diabetes, link, pathfinder, pigs, munin1, munin2, munin3, munin4';

% Massive Networks (>1000 nodes)   munin

% ----- Continues Bayesian Networks (Gaussian)

% Medium Networks (20每50 nodes) ecoli70, magic-niab

% Large Networks (50每100 nodes) magic-irri

% Very Large Networks (101每1000 nodes) arth150

% ----- Continues Bayesian Networks (Conditional Linear Gaussian)

% Small Networks (<20 nodes) sangiovese

% Medium Networks (20每50 nodes) mehra-original, mehra-complete


% ************ Names of the 7 algorithms are as follows ************ 

% ----- Global causal structure learning  GSBN, GES, PC, MMHC, PCstable, F2SL_c, F2SL_s



clear all
clc
close all


% Name of data
data_name='alarm';

% Samples of data
data_samples=5000;

% 'dis' represents discrete data, 'con' denotes continues data 
data_type='dis';      

% Name of algorithm
alg_name='MMHC';

% Significance level
alpha=0.01;


% Path of the data set
data_path=strcat('data/',data_name,'_',num2str(data_samples),'.txt');
if exist(data_path,'file')==0
     fprintf('\n%s does not exist.\n\n',strcat(data_name,'_',num2str(data_samples),'.txt'));
     return;
end

% Load data according to the path
% data needs to start from 0
data = importdata(data_path)+1;

% data(:,16)=[];


% Causal_Learner
[Result1,Result2]=Causal_Learner(alg_name,data,data_type,alpha);

% Global causal structure learning 
% Result1 is learned causal graph
% Result2 is running time



% -------------------------------------
% -------------------------------------
% Evaluation


% Path of the graph
graph_path=strcat('data/',data_name,'_graph.txt');
if exist(graph_path,'file')==0
     fprintf('\n%s does not exist.\n\n',strcat(data_name,'_graph.txt'));
     return;
end

% Load graph (true DAG) according to the path
graph = importdata(graph_path);


% -------------------------------------
% Evaluate global causal structure

DAG=Result1;
[arrhd_F1,arrhd_precision,arrhd_recall,SHD,reverse,miss,extra]=evaluation_GCS(DAG,graph);
fprintf('\nThe learned global causal structure is as follows.\n');
sparse(DAG)
draw_graph(DAG);
fprintf('arrhd_F1=%.2f, arrhd_precision=%.2f, arrhd_recall=%.2f\n\n',arrhd_F1,arrhd_precision,arrhd_recall);
fprintf('SHD=%.0f, reverse=%.0f, miss=%.0f, extra=%.0f\n',SHD,reverse,miss,extra);
fprintf('\nElapsed time is %.2f seconds.\n\n\n',Result2);





