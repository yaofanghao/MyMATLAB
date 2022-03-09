%--------------------------------------------------------------------------
%           创建一个Microsoft Excel服务器，改变单元格A1的边框颜色
%--------------------------------------------------------------------------
% CopyRight：xiezhh

Excel = actxserver('Excel.Application');
Excel.Visible = 1;
Workbook = Excel.Workbooks.Add;
Sheet1 = Workbook.Sheets.Item(1);

% 通过循环改变单元格A1的边框颜色
for i = 0:56
    Sheet1.Range('A1').Borders.ColorIndex = i;
    pause(1);
end