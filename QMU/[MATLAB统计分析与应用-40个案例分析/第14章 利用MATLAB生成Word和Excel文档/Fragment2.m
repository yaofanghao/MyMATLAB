%--------------------------------------------------------------------------
%               ����һ��Microsoft Excel������������ͼ14.4-1
%--------------------------------------------------------------------------
% CopyRight��xiezhh

Excel = actxserver('Excel.Application');
Excel.Visible = 1;
Workbook = Excel.Workbooks.Add;
Sheet1 = Workbook.Sheets.Item(1);

for j=1:2:7
    for i=1:14
        Sheet1.Cells.Item((i-1)*256 + j).Interior.ColorIndex = i + 14*(j-1)/2;
        Sheet1.Cells.Item((i-1)*256 + j + 1).Value = i + 14*(j-1)/2;
        Sheet1.Cells.Item((i-1)*256 + j + 1).HorizontalAlignment = 2;
    end
end

% ע��������������Excel2003�汾������Excel2007��Excel2010�����У���Ҫ��256��Ϊ16384