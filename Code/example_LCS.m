

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


% ************ Names of the 4 algorithms are as follows ************ 

% ----- Local causal structure learning   PCD_by_PCD, MB_by_MB, CMB, LCS_FS


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
alg_name='ELCS';

% Significance level
alpha=0.01;

% Index of target node. If it is global structure learning, this parameter is not needed
target=23;   

% Path of the data set
data_path=strcat('data/',data_name,'_',num2str(data_samples),'.txt');
if exist(data_path,'file')==0
     fprintf('\n%s does not exist.\n\n',strcat(data_name,'_',num2str(data_samples),'.txt'));
     return;
end

% Load data according to the path
% data needs to start from 0
data = importdata(data_path)+1;


% Causal_Learner
[Result1,Result2,Result3,Result4,Result5]=Causal_Learner(alg_name,data,data_type,alpha,target);

% Local causal structure learning
% Result1 is learned target's parents.
% Result2 is learned target's children.
% Result3 is learned target's PC, but cannot distinguish whether they are parents or children.
% Result4 is the number of conditional independence tests
% Result5 is running time


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
% Evaluate local causal structure

Parents=Result1;
Children=Result2;
Undirected=Result3;
[arrhd_F1,arrhd_precision,arrhd_recall,SHD,reverse,miss,extra,undirected]=evaluation_LCS(Parents,Children,Undirected,target,graph);
fprintf('\nThe learned parents of target %.0f is [',target);
for i=1:length(Parents)
    if i==length(Parents)
        fprintf('%d',Parents(i));
    else
        fprintf('%d\t',Parents(i));
    end
end
fprintf(']\n\nThe learned children of target %.0f is [',target);
for i=1:length(Children)
    if i==length(Children)
        fprintf('%d',Children(i));
    else
        fprintf('%d\t',Children(i));
    end
end
fprintf(']\n\narrhd_F1=%.2f, arrhd_precision=%.2f, arrhd_recall=%.2f\n\n',arrhd_F1,arrhd_precision,arrhd_recall);
fprintf('SHD=%.0f, reverse=%.0f, miss=%.0f, extra=%.0f, undirected=%.0f\n',SHD,reverse,miss,extra,undirected);
fprintf('\nThe number of conditional independence tests is %.0f.\n',Result4);
fprintf('\nElapsed time is %.2f seconds.\n\n\n',Result5);
    










