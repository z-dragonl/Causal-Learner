

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


% ************ Names of the 15 algorithms are as follows ************ 

% ----- Markov blanket learning   GS, IAMB, interIAMB, IAMBnPC, interIAMBnPC, FastIAMB, FBED
%                                 MMMB, HITONMB, PCMB, IPCMB, MBOR, STMB, BAMB, EEMB, MBFS


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
alg_name='HITONMB';

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
[Result1,Result2,Result3]=Causal_Learner(alg_name,data,data_type,alpha,target);

% Markov blanket learning 
% Result1 is learned target's Markov blanket.
% Result2 is the number of conditional independence tests
% Result3 is running time



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
% Evaluate Markov blanket

MB=Result1;
[adj_F1,adj_precision,adj_recall]=evaluation_MB(MB,target,graph);
fprintf('\nThe learned Markov blanket of target %.0f is [',target);
for i=1:length(MB)
    if i==length(MB)
        fprintf('%d',MB(i));
    else
        fprintf('%d\t',MB(i));
    end
end
fprintf(']\n\nadj_F1=%.2f, adj_precision=%.2f, adj_recall=%.2f\n',adj_F1,adj_precision,adj_recall);
fprintf('\nThe number of conditional independence tests is %.0f.\n',Result2);
fprintf('\nElapsed time is %.2f seconds.\n\n\n',Result3);
















