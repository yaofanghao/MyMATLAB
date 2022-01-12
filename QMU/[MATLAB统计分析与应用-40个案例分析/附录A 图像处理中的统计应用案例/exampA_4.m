%--------------------------------------------------------------------------
%  exampA.4-1  基于贝叶斯判别的手写体数字识别
%--------------------------------------------------------------------------
% CopyRight：xiezhh

%***********************读取训练样本图片，并进行预处理***********************
% 利用uigetfile函数交互式选取训练样本图片
geshi = {'*.jpg','JPEG image (*.jpg)';...
       '*.bmp','Bitmap image (*.bmp)';...
       '*.*','All Files (*.*)'};
[FileName FilePath] = uigetfile(geshi,'选取训练样本图片',...
    '*.jpg','MultiSelect','on');

if ~isequal([FileName,FilePath],[0,0]);
    FileFullName = strcat(FilePath,FileName);
    if ~iscell(FileName)
        FileName = {FileName};
        FileFullName = {FileFullName};
    end
else
    return;
end
n = length(FileFullName);
I = zeros(50);
BW = zeros(16);
training = zeros(1,256);
group = [];
% 通过循环对每一个图片进行标准化处理，并生成训练样本数据矩阵training和分组向量group
for i = 1:n
    I = imread(FileFullName{i});
    I = 255-I;   
    I = im2bw(I,0.4);
    [y,x] = find(I == 1);
    BW = I(min(y):max(y),min(x):max(x));
    BW = imresize(BW,[16, 16]);
    training(i,:) = double(BW(:)');
    [pathstr,namestr,ext] = fileparts(FileName{i});
    group = [group;str2num(namestr(4))];
end


%***************************创建朴素贝叶斯分类器对象*************************
% 创建朴素贝叶斯分类器对象ObjBayes
ObjBayes = NaiveBayes.fit(training,group,'Distribution','mn')


%*****************************对训练样本图片进行判别*************************
% 利用所创建的朴素贝叶斯分类器对象ObjBayes，对训练样本图片进行判别
pre0 = ObjBayes.predict(training);
isequal(pre0, group)


%***********************读取检验样本图片，并进行预处理***********************
% 利用uigetfile函数交互式选取检验样本图片
geshi = {'*.jpg','JPEG image (*.jpg)';...
       '*.bmp','Bitmap image (*.bmp)';...
       '*.*','All Files (*.*)'};
[FileName FilePath] = uigetfile(geshi,'选取检验样本图片',...
    '*.jpg','MultiSelect','on');

if ~isequal([FileName,FilePath],[0,0]);
    FileFullName = strcat(FilePath,FileName);
    if ~iscell(FileName)
        FileName = {FileName};
        FileFullName = {FileFullName};
    end
else
    return;
end
n = length(FileFullName);
I = zeros(50);
BW = zeros(16);
sampledata = zeros(1,256);
samplegroup = [];
% 通过循环对每一个图片进行标准化处理，并生成检验样本数据矩阵sampledata和分组向量samplegroup
for i = 1:n
    I = imread(FileFullName{i});
    I = 255-I;
    I = im2bw(I,0.4);
    [y,x] = find(I == 1);
    BW = I(min(y):max(y),min(x):max(x));
    BW = imresize(BW,[16, 16]);
    sampledata(i,:) = double(BW(:)');
    [pathstr,namestr,ext] = fileparts(FileName{i});
    samplegroup = [samplegroup; str2num(namestr(4))];
end


%*****************************对检验样本图片进行判别*************************
% 利用所创建的朴素贝叶斯分类器对象ObjBayes，对检验样本图片进行判别
pre1 = ObjBayes.predict(sampledata);
% 查看判别结果
[samplegroup, pre1]