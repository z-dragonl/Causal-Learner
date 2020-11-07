
function [Result1,Result2,Result3,Result4,Result5] = Causal_Learner(input_alg_name,data,data_type,alpha,target)

if (nargin < 4)
    fprintf('\nInput parameters is not valid.\n\n');
    return;
end

addpath(genpath(pwd));

% Maximum size of conditioning set
maxK=3;

% Get the number of samples in the data set: samples
% and the number of nodes: p
[samples,p]=size(data);

% Size array of each node
ns=max(data);

total_alg_num=27;

total_alg=cell(1,total_alg_num);

for alg_i=1:total_alg_num
    
    % Global causal structure learning
    
    if alg_i==1
        alg_name='GSBN';
    elseif alg_i==2
        alg_name='GES';
    elseif alg_i==3
        alg_name='PC';
    elseif alg_i==4
        alg_name='MMHC';
    elseif alg_i==5
        alg_name='PCstable';
    elseif alg_i==6
        alg_name='F2SL_c';
    elseif alg_i==7
        alg_name='F2SL_s';      
        
    % Local causal structure learning
        
    elseif alg_i==8
        alg_name='PCD_by_PCD';
    elseif alg_i==9
        alg_name='MB_by_MB';
    elseif alg_i==10
        alg_name='CMB';
    elseif alg_i==11
        alg_name='LCS_FS';
        
    % Markov blanket learning       
    
    elseif alg_i==12
        alg_name='GS';
    elseif alg_i==13
        alg_name='IAMB';
    elseif alg_i==14
        alg_name='interIAMB';
    elseif alg_i==15
        alg_name='IAMBnPC';
    elseif alg_i==16
        alg_name='interIAMBnPC';
    elseif alg_i==17
        alg_name='FastIAMB';
    elseif alg_i==28
        alg_name='FBED';
        
    elseif alg_i==19
        alg_name='MMMB';
    elseif alg_i==20
        alg_name='HITONMB';
    elseif alg_i==21
        alg_name='PCMB';
    elseif alg_i==22
        alg_name='IPCMB';
    elseif alg_i==23
        alg_name='MBOR';
    elseif alg_i==24
        alg_name='STMB';
    elseif alg_i==25
        alg_name='BAMB';
    elseif alg_i==26
        alg_name='EEMB';
    elseif alg_i==27
        alg_name='MBFS';
    end
    
    total_alg{alg_i}=alg_name;
    
end

    
if ~ismember(input_alg_name,total_alg)
    fprintf('\n%s is not a valid algorithm name \n',input_alg_name);
    return;
end


% Identify the type of algorithm

alg_index=strmatch(input_alg_name,total_alg);

isGCS=0;
isLCS=0;
isMB=0;

if alg_index<=7
    isGCS=1;
elseif alg_index<=11
    isLCS=1;
elseif alg_index<=27
    isMB=1;
end


if isGCS
    fprintf('\nGlobal causal structure learning by %s\n\n',input_alg_name);
elseif isLCS
    fprintf('\nLocal causal structure learning by %s\n\n',input_alg_name);
elseif isMB
    fprintf('\nMarkov blanket learning by %s\n\n',input_alg_name);
end


if strcmp(data_type,'dis')
    
    % using G2 test
    algorithm=str2func(strcat(input_alg_name,'_G2'));
    
    if isGCS
        [DAG,time]=algorithm(data,alpha,ns,p,maxK);
    elseif isLCS
        [Parents,Children,PC,Undirected,test,time] = algorithm(data,target,alpha,ns,p,maxK);
    elseif isMB
        [MB,test,time] = algorithm (data,target,alpha,ns,p,maxK);
    end
        
elseif strcmp(data_type,'con')
    
    if strcmp(input_alg_name,'FastIAMB')||strcmp(input_alg_name,'MMHC')||strcmp(input_alg_name,'F2SL_s')
        fprintf('\nSorry, %s only supports discrete data.\n\n',input_alg_name);
        return;
    end
    
    % using Fisher's Z test
    algorithm=str2func(strcat(input_alg_name,'_Z'));
    
    if isGCS
        [DAG,time]=algorithm(data,alpha,samples,p,maxK);
    elseif isLCS
        [Parents,Children,PC,Undirected,test,time] = algorithm(data,target,alpha,samples,p,maxK);
    elseif isMB
        [MB,test,time] = algorithm (data,target,alpha,samples,p,maxK);
    end
end

    
fprintf('\n------------------------------------\n\n');

if isGCS    
    Result1=DAG;
    Result2=time;
    Result3=[];
    Result4=[];
    Result5=[];
elseif isLCS  
    Result1=Parents;
    Result2=Children;
    Result3=Undirected;
    Result4=test;
    Result5=time;
elseif isMB  
    Result1=MB;
    Result2=test;
    Result3=time;
    Result4=[];
    Result5=[];
end
    
    
    
    
    
    
    
    
    
    

