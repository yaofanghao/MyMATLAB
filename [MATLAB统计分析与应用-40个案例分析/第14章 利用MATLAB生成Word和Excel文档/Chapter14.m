%--------------------------------------------------------------------------
%  第14章 利用MATLAB生成Word和Excel文档
%--------------------------------------------------------------------------
% CopyRight：xiezhh

%%  examp14.2-1  创建一个日历控件
f = figure('position', [360  278  535  410]); 
cal = actxcontrol('mscal.calendar', [0 0 535 410], f) 

%%  examp14.2-2  invoke、interfaces 、methods、methodsview和events函数的用法示例
f = figure('position', [100 100 600 500]);
cal = actxcontrol('mscal.calendar', [0 0 600 500], f);
cal.invoke
cal.interfaces
cal.methods
cal.events
cal.methodsview

%%  examp14.2-3  iscom、isinterface、isprop、ismethod和isevent函数的用法示例
Word = actxserver('Word.Application');
Word.iscom
Word.Documents.isinterface
isprop(Word, 'Width')
ismethod(Word, 'Quit')
isevent(Word, 'Quit')

%%  examp14.2-4  get和inspect函数的用法示例
Word = actxserver('Word.Application');
get(Word)
Word.get
Word.Visible
Word.get('Visible')
get(Word, 'Visible')
Word.inspect

%%  examp14.2-5  属性传参示例
Word = actxserver('Word.Application'); 
set(Word, 'Visible', 1);
document = invoke(Word.Documents, 'Add'); 
Content = document.Content;
isinterface(Content)
end_of_doc = get(Content,'end')

%%  examp14.2-6  set、addproperty和deleteproperty函数的用法示例
f = figure('position', [100 100 600 500]);
cal = actxcontrol('mscal.calendar', [0 0 600 500], f);
cal.set('month',10)

h = actxcontrol('mwsamp.mwsampctrl.2',[200 120 200 200]);
h.addproperty('xiezhh');
h.xiezhh = 'I''m a teacher';
h.get
h.deleteproperty('xiezhh');
h.get

%%  examp14.2-7  创建一个日历控件，调用NextDay方法修改时间
figure;
cal = actxcontrol('mscal.calendar',[10 10 540 400]);
for i = 1:1000
       cal.NextDay;
end
cal.Value

%%  examp14.2-8  创建一个日历控件，为控件注册通用事件处理程序
figure;
cal = actxcontrol('mscal.calendar',[10 10 540 400]);
cal.events
cal.registerevent('XiezhhTest');

%%  examp14.2-9  创建一个日历控件，保存其初始状态，稍后重新加载
figure;
cal = actxcontrol('mscal.calendar',[10 10 540 400]);
cal.save('mscal.mat');

cal.month = 1;
cal.day = 1;
cal.year = 2000;
cal.get

cal.load('mscal.mat');    % 重新加载cal的初始状态
cal.get

%%  examp14.2-10  创建Windows Media Player服务器，播放歌曲“青花瓷”
h = actxserver('WMPlayer.OCX.7')
h.get
h.invoke
h.openPlayer('F:\我的音乐盒\青花瓷.mp3')

%%  examp14.2-11  创建一个运行Microsoft  Excel的COM服务器应用程序
Excel = actxserver('Excel.Application'); 
set(Excel, 'Visible', 1);
Workbooks = Excel.Workbooks;
Workbook = Workbooks.invoke('Add');

%%  examp14.3-1 和 examp14.3-2  创建一个Microsoft Word服务器，插入图片
% 创建一个Microsoft Word服务器
Word = actxserver('Word.Application'); 
Word.Visible = 1;
Document = Word.Documents.Add;
Selection = Word. Selection;

%********************************插入外部图片*******************************
filename = [matlabroot '\toolbox\images\imdemos\football.jpg'];
handle1 = Selection.InlineShapes.AddPicture(filename);
handle2 = Document.Shapes.AddPicture(filename, [], [], 180, 50, 200, 170);

%********************************插入内部图片*******************************
Document = Word.Documents.Add;
Selection = Word. Selection;

rng('default');
data = normrnd(75,6,1000,1);

zft = figure('units','normalized','position',...
                [0.280469 0.553385 0.428906 0.251302],'visible','off');
set(gca,'position',[0.1 0.2 0.85 0.75]);
hist(data);
grid on;
xlabel('考试成绩');
ylabel('人数');
hgexport(zft, '-clipboard');
delete(zft);

Selection.Paste;
Selection.TypeParagraph;
Selection.PasteSpecial;