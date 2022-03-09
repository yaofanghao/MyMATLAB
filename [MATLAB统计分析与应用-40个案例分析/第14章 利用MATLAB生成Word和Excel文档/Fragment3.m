%--------------------------------------------------------------------------
%               创建一个Microsoft Excel服务器，生成图14.4-2
%--------------------------------------------------------------------------
% CopyRight：xiezhh

Excel = actxserver('Excel.Application');
Excel.Visible = 1;
Workbook = Excel.Workbooks.Add;
Sheet1 = Workbook.Sheets.Item(1);
Tabhead = {'样式','LineStyle','Weight','样式','LineStyle','Weight'};
Sheet1.Range('A1:F1').Value = Tabhead;
Sheet1.Range('A2').Value = '无';

for j = 1:3:4
    for i=2:8
        Sheet1.Cells.Item((i-1)*256+j).Borders.Item(4).Linestyle = i-2+7*(j-1)/3;
        Weight = Sheet1.Cells.Item((i-1)*256 + j).Borders.Item(4).Weight;
        Sheet1.Cells.Item((i-1)*256+j).Borders.Item(4).ColorIndex = 1;
        Sheet1.Cells.Item((i-1)*256+j+1).Value = i-2+7*(j-1)/3;
        Sheet1.Cells.Item((i-1)*256+j+1).HorizontalAlignment = 2;     
        Sheet1.Cells.Item((i-1)*256+j+2).Value = Weight;
        Sheet1.Cells.Item((i-1)*256+j+2).HorizontalAlignment = 2;
    end
end

% 注：本程序适用于Excel2003版本，若在Excel2007或Excel2010下运行，需要将256改为16384