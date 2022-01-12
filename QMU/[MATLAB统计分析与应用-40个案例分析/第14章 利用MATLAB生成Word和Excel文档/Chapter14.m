%--------------------------------------------------------------------------
%  ��14�� ����MATLAB����Word��Excel�ĵ�
%--------------------------------------------------------------------------
% CopyRight��xiezhh

%%  examp14.2-1  ����һ�������ؼ�
f = figure('position', [360  278  535  410]); 
cal = actxcontrol('mscal.calendar', [0 0 535 410], f) 

%%  examp14.2-2  invoke��interfaces ��methods��methodsview��events�������÷�ʾ��
f = figure('position', [100 100 600 500]);
cal = actxcontrol('mscal.calendar', [0 0 600 500], f);
cal.invoke
cal.interfaces
cal.methods
cal.events
cal.methodsview

%%  examp14.2-3  iscom��isinterface��isprop��ismethod��isevent�������÷�ʾ��
Word = actxserver('Word.Application');
Word.iscom
Word.Documents.isinterface
isprop(Word, 'Width')
ismethod(Word, 'Quit')
isevent(Word, 'Quit')

%%  examp14.2-4  get��inspect�������÷�ʾ��
Word = actxserver('Word.Application');
get(Word)
Word.get
Word.Visible
Word.get('Visible')
get(Word, 'Visible')
Word.inspect

%%  examp14.2-5  ���Դ���ʾ��
Word = actxserver('Word.Application'); 
set(Word, 'Visible', 1);
document = invoke(Word.Documents, 'Add'); 
Content = document.Content;
isinterface(Content)
end_of_doc = get(Content,'end')

%%  examp14.2-6  set��addproperty��deleteproperty�������÷�ʾ��
f = figure('position', [100 100 600 500]);
cal = actxcontrol('mscal.calendar', [0 0 600 500], f);
cal.set('month',10)

h = actxcontrol('mwsamp.mwsampctrl.2',[200 120 200 200]);
h.addproperty('xiezhh');
h.xiezhh = 'I''m a teacher';
h.get
h.deleteproperty('xiezhh');
h.get

%%  examp14.2-7  ����һ�������ؼ�������NextDay�����޸�ʱ��
figure;
cal = actxcontrol('mscal.calendar',[10 10 540 400]);
for i = 1:1000
       cal.NextDay;
end
cal.Value

%%  examp14.2-8  ����һ�������ؼ���Ϊ�ؼ�ע��ͨ���¼��������
figure;
cal = actxcontrol('mscal.calendar',[10 10 540 400]);
cal.events
cal.registerevent('XiezhhTest');

%%  examp14.2-9  ����һ�������ؼ����������ʼ״̬���Ժ����¼���
figure;
cal = actxcontrol('mscal.calendar',[10 10 540 400]);
cal.save('mscal.mat');

cal.month = 1;
cal.day = 1;
cal.year = 2000;
cal.get

cal.load('mscal.mat');    % ���¼���cal�ĳ�ʼ״̬
cal.get

%%  examp14.2-10  ����Windows Media Player�����������Ÿ������໨�ɡ�
h = actxserver('WMPlayer.OCX.7')
h.get
h.invoke
h.openPlayer('F:\�ҵ����ֺ�\�໨��.mp3')

%%  examp14.2-11  ����һ������Microsoft  Excel��COM������Ӧ�ó���
Excel = actxserver('Excel.Application'); 
set(Excel, 'Visible', 1);
Workbooks = Excel.Workbooks;
Workbook = Workbooks.invoke('Add');

%%  examp14.3-1 �� examp14.3-2  ����һ��Microsoft Word������������ͼƬ
% ����һ��Microsoft Word������
Word = actxserver('Word.Application'); 
Word.Visible = 1;
Document = Word.Documents.Add;
Selection = Word. Selection;

%********************************�����ⲿͼƬ*******************************
filename = [matlabroot '\toolbox\images\imdemos\football.jpg'];
handle1 = Selection.InlineShapes.AddPicture(filename);
handle2 = Document.Shapes.AddPicture(filename, [], [], 180, 50, 200, 170);

%********************************�����ڲ�ͼƬ*******************************
Document = Word.Documents.Add;
Selection = Word. Selection;

rng('default');
data = normrnd(75,6,1000,1);

zft = figure('units','normalized','position',...
                [0.280469 0.553385 0.428906 0.251302],'visible','off');
set(gca,'position',[0.1 0.2 0.85 0.75]);
hist(data);
grid on;
xlabel('���Գɼ�');
ylabel('����');
hgexport(zft, '-clipboard');
delete(zft);

Selection.Paste;
Selection.TypeParagraph;
Selection.PasteSpecial;