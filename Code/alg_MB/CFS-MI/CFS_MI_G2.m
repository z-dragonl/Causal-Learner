
function   [mb,ntest,time]=CFS_MI_G2(Data,target,alpha,ns,p,k)
start=tic;
ntest=0;
mb=[];
% time=0;
% 求PC
% pc = yingshe2(pc,target)
%ntest=ntest+ntest1;
train_data=Data(:,mysetdiff(1:p,target));
featureMatrix = train_data;
% 标签，这里就是target这列的数据
train_label=Data(:,target);
classColumn = train_label;
numFeatures = size(featureMatrix,2);
% classScore构造一个全零的矩阵有numFeatures行，一列
classScore = zeros(numFeatures,1);
vis = zeros(p,1);
for i = 1:numFeatures
    ntest=ntest+1;
%     iXY(i) = mi(featureMatrix(:,i),classColumn);
    classScore(i) = SU(featureMatrix(:,i),classColumn);
end
% 排序，1--列的意思，desend降序，indexScore排序后元素在原矩阵中的位置
[classScore, indexScore] = sort(classScore,1,'descend');
% [iXY, iXYzhi] = sort(iXY,1,'descend');

% th2 = 0.0001;
% 这个值用来控制进去CPC的数量，越小进去的越多
threshold = 0.05;


th3 = 0.15;
t = indexScore(classScore < threshold);
u = classScore(classScore < threshold);
% 更新两个矩阵，保留classScore>threhold的元素
indexScore = indexScore(classScore > threshold);
classScore = classScore(classScore > threshold);
% 判断更新后的矩阵是否为空，为空直接结束函数，选择的PC也就是空，不空的话设置一个当前位置变量为1，接下来用这个变量作为循环控制变量
if ~isempty(indexScore)
    curPosition = 1;
else
    curPosition = 0;
    selectedFeatures=[];
    return;
end
% 循环开始   

mii = -1;
while curPosition <= length(indexScore)
    mb_tmp = [];
    j = curPosition + 1;
    % curFeature当前排在最前面的一个点
    curFeature = indexScore(curPosition);
    while j <= length(indexScore)
        
        % 计算XY之间的SU
        scoreij = SU(featureMatrix(:,curFeature),featureMatrix(:,indexScore(j)));
        ntest=ntest+1;
        % 如果XY之间的SU比Y和T之间的SU要大，Y就可以删除了，这个证明了的
        if scoreij > classScore(j)
            indexScore(j) = [];
            classScore(j) = [];
%             mb_tmp = myunion(mb_tmp,indexScore(j));
        else
%             mii = classScore(j);
            j = j + 1;
        end
    end
    curPosition = curPosition + 1;
end
selectedFeatures = indexScore;
pc = selectedFeatures;
last = selectedFeatures(end);
mb = yingshe2(pc,target);
for i =1:length(selectedFeatures)
    vis(selectedFeatures(i))=1;
end
% 配偶



len1 = length(selectedFeatures);
for i=1:len1
    mb_tmp=[];
    a = find(selectedFeatures == last)+1;
    len2 = length(t);
    while a <= len2
        
        if vis(t(a))==1%ismember(t(a),mb_tmp)
            a = a + 1;
            continue;
        end
%         if ~vis(selectedFeatures(i),t(a))
%             vis(selectedFeatures(i),t(a)) = SU(featureMatrix(:,selectedFeatures(i)),featureMatrix(:,t(a)));
%             vis(t(a),selectedFeatures(i)) = vis(selectedFeatures(i),t(a));
%         end
        scoreij = SU(featureMatrix(:,selectedFeatures(i)),featureMatrix(:,t(a)));   
        ntest=ntest+1;
        % 这个判断用来控制真正进去CPC的数目，阈值越低进去的越多，但是进去的多可能有错误节点
        if scoreij > u(a)+ 0.13
%             mb_tmp = myunion(mb_tmp,t(a));
            iXYZ = cmi(train_label,featureMatrix(:,t(a)),featureMatrix(:,selectedFeatures(i)));
            iXY0 = mi(train_label,featureMatrix(:,t(a)));
            ntest=ntest+2;
            if t(a)>=target
                ttt = t(a)+1;
            else
                ttt = t(a);
            end 
            if iXYZ > iXY0
                mb_tmp = myunion(mb_tmp,t(a));
                mb = myunion(mb,ttt);
                vis(t(a))=1;
            end    
        end
        a = a + 1;
    end
end


time=toc(start);

end

function [score] = SU(firstVector,secondVector)
%function [score] = SU(firstVector,secondVector)
%calculates SU = 2 * (I(X;Y)/(H(X) + H(Y)))
hX = h(firstVector);
hY = h(secondVector);
iXY = mi(firstVector,secondVector);

score = (2 * iXY) / (hX + hY);
end


function pc = yingshe2( pc ,target)
%YINGSHE2 此处显示有关此函数的摘要
%   此处显示详细说明
for i=1:length(pc)
    if pc(i)>=target
        pc(i)= pc(i)+1;
    end
end

end

        






