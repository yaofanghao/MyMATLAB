function [newData1] = importfile(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  从指定文件中导入数据
%  FILETOREAD1:  要读取的文件

%  由 MATLAB 于 29-May-2022 23:49:33 自动生成

% 导入文件
sheetName='Sheet1';
[numbers, strings, raw] = xlsread(fileToRead1, sheetName);
if ~isempty(numbers)
    newData1.data =  numbers;
end
if ~isempty(strings)
    newData1.textdata =  strings;
end

if ~isempty(strings) && ~isempty(numbers)
    [strRows, strCols] = size(strings);
    [numRows, numCols] = size(numbers);
    likelyRow = size(raw,1) - numRows;
    % 将数据拆分为每列包含一个字段的新结构体。
    if strCols == numCols && likelyRow > 0 && strRows >= likelyRow
        newData1.colheaders = strings(likelyRow, :);
    end
end

