%--------------------------------------------------------------------------
%           ����һ��Microsoft Excel���������ı䵥Ԫ��A1�ı߿���ɫ
%--------------------------------------------------------------------------
% CopyRight��xiezhh

Excel = actxserver('Excel.Application');
Excel.Visible = 1;
Workbook = Excel.Workbooks.Add;
Sheet1 = Workbook.Sheets.Item(1);

% ͨ��ѭ���ı䵥Ԫ��A1�ı߿���ɫ
for i = 0:56
    Sheet1.Range('A1').Borders.ColorIndex = i;
    pause(1);
end