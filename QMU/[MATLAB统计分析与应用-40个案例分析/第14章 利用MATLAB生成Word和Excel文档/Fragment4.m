%--------------------------------------------------------------------------
%               创建一个Microsoft Excel 服务器，插入图片
%--------------------------------------------------------------------------
% CopyRight：xiezhh

% 创建一个Microsoft Excel 服务器
Excel = actxserver('Excel.Application');
Excel.Visible = 1;
Workbook = Excel.Workbooks.Add;


%********************************插入外部图片*******************************
PicturePath = which('peppers.png');
% 在当前工作表的指定位置处插入一幅指定大小的图片
h1 = Excel.ActiveSheet.Shapes.AddPicture(PicturePath,0,1,50,60,400,300);
h2 = Excel.ActiveSheet.Shapes.AddPicture(PicturePath,0,1,500,60,200,150);
h3 = Excel.ActiveSheet.Shapes.AddPicture(PicturePath,0,1,650,180,200,150);


%********************************插入内部图片*******************************
Workbook = Excel.Workbooks.Add;
Sheet1 = Workbook.Sheets.Item(1);

% 生成随机数，绘制直方图
data = normrnd(75,4,1000,1);
zft = figure('units','normalized','position',...
                [0.280469 0.553385 0.428906 0.251302],'visible','off');
set(gca,'position',[0.1 0.2 0.85 0.75]);
hist(data);
grid on;
xlabel('考试成绩');
ylabel('人数');
hgexport(zft, '-clipboard');

% 选中工作表Sheet1的A2单元格，插入由MATLAB命令作出的直方图
Sheet1.Range('A2').Select;
Sheet1.Paste

Sheet1.Range('E11').Select;
Sheet1.PasteSpecial;

Sheet1. Range('I20').PasteSpecial;