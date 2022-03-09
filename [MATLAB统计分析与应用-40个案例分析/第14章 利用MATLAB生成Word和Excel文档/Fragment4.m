%--------------------------------------------------------------------------
%               ����һ��Microsoft Excel ������������ͼƬ
%--------------------------------------------------------------------------
% CopyRight��xiezhh

% ����һ��Microsoft Excel ������
Excel = actxserver('Excel.Application');
Excel.Visible = 1;
Workbook = Excel.Workbooks.Add;


%********************************�����ⲿͼƬ*******************************
PicturePath = which('peppers.png');
% �ڵ�ǰ�������ָ��λ�ô�����һ��ָ����С��ͼƬ
h1 = Excel.ActiveSheet.Shapes.AddPicture(PicturePath,0,1,50,60,400,300);
h2 = Excel.ActiveSheet.Shapes.AddPicture(PicturePath,0,1,500,60,200,150);
h3 = Excel.ActiveSheet.Shapes.AddPicture(PicturePath,0,1,650,180,200,150);


%********************************�����ڲ�ͼƬ*******************************
Workbook = Excel.Workbooks.Add;
Sheet1 = Workbook.Sheets.Item(1);

% ���������������ֱ��ͼ
data = normrnd(75,4,1000,1);
zft = figure('units','normalized','position',...
                [0.280469 0.553385 0.428906 0.251302],'visible','off');
set(gca,'position',[0.1 0.2 0.85 0.75]);
hist(data);
grid on;
xlabel('���Գɼ�');
ylabel('����');
hgexport(zft, '-clipboard');

% ѡ�й�����Sheet1��A2��Ԫ�񣬲�����MATLAB����������ֱ��ͼ
Sheet1.Range('A2').Select;
Sheet1.Paste

Sheet1.Range('E11').Select;
Sheet1.PasteSpecial;

Sheet1. Range('I20').PasteSpecial;