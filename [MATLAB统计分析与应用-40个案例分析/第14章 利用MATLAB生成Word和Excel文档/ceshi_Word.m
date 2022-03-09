function ceshi_Word
% 利用MATLAB生成Word文档
%
% CopyRight：xiezhh

filespec_user = [pwd '\测试.doc'];

% 判断Word是否已经打开，若已打开，就在打开的Word中进行操作，否则就打开Word
try
    Word = actxGetRunningServer('Word.Application');
catch
    Word = actxserver('Word.Application'); 
end;

Word.Visible = 1;

% 若测试文件存在，打开该测试文件，否则，新建一个文件，并保存，文件名为测试.doc
if exist(filespec_user,'file'); 
    Document = Word.Documents.Open(filespec_user);     
    % Document = invoke(Word.Documents,'Open',filespec_user);
else
    Document = Word.Documents.Add;
    try
        Document.SaveAs(filespec_user);
    catch
        Document.SaveAs2(filespec_user);
    end
end

Content = Document.Content;
Selection = Word.Selection;
Paragraphformat = Selection.ParagraphFormat;

% 页面设置
Document.PageSetup.TopMargin = 60;
Document.PageSetup.BottomMargin = 45;
Document.PageSetup.LeftMargin = 45;
Document.PageSetup.RightMargin = 45;

% 设定文档内容的起始位置和标题
Content.Start = 0;
headline = '试  卷  分  析';
Content.Text = headline;
Content.Font.Size = 16 ;
Content.Font.Bold = 4 ;
Content.Paragraphs.Alignment = 'wdAlignParagraphCenter';

Selection.Start = Content.end;
Selection.TypeParagraph;

xueqi = '（ 2009  —  2010   学年 第一学期）';
Selection.Text = xueqi;
Selection.Font.Size = 12;
Selection.Font.Bold = 0;
Selection.MoveDown;
Paragraphformat.Alignment = 'wdAlignParagraphCenter';
Selection.TypeParagraph;
Selection.TypeParagraph;
Selection.Font.Size = 10.5;

% 在光标所在位置插入一个12行9列的表格
Tables = Document.Tables.Add(Selection.Range,12,9);    

DTI = Document.Tables.Item(1);          % 或DTI = Tables;

% 设置表格边框
DTI.Borders.OutsideLineStyle = 'wdLineStyleSingle';
DTI.Borders.OutsideLineWidth = 'wdLineWidth150pt';
DTI.Borders.InsideLineStyle = 'wdLineStyleSingle';
DTI.Borders.InsideLineWidth = 'wdLineWidth150pt';
DTI.Rows.Alignment = 'wdAlignRowCenter';
DTI.Rows.Item(8).Borders.Item(1).LineStyle = 'wdLineStyleNone';
DTI.Rows.Item(8).Borders.Item(3).LineStyle = 'wdLineStyleNone';
DTI.Rows.Item(11).Borders.Item(1).LineStyle = 'wdLineStyleNone';
DTI.Rows.Item(11).Borders.Item(3).LineStyle = 'wdLineStyleNone';

% 设置表格列宽和行高
column_width = [53.7736,85.1434,53.7736,35.0094,...
    35.0094,76.6981,55.1887,52.9245,54.9057];
row_height = [28.5849,28.5849,28.5849,28.5849,25.4717,25.4717,...
    32.8302,312.1698,17.8302,49.2453,14.1509,18.6792];
for i = 1:9
    DTI.Columns.Item(i).Width = column_width(i);
end
for i = 1:12
    DTI.Rows.Item(i).Height = row_height(i);
end

% 通过循环设置每个单元格的垂直对齐方式
for i = 1:12
    for j = 1:9
        DTI.Cell(i,j).VerticalAlignment = 'wdCellAlignVerticalCenter';
    end
end

% 合并单元格
DTI.Cell(1, 4).Merge(DTI.Cell(1, 5));
DTI.Cell(2, 4).Merge(DTI.Cell(2, 5));
DTI.Cell(3, 4).Merge(DTI.Cell(3, 5));
DTI.Cell(4, 4).Merge(DTI.Cell(4, 5));
DTI.Cell(5, 2).Merge(DTI.Cell(5, 5));
DTI.Cell(5, 3).Merge(DTI.Cell(5, 6));
DTI.Cell(6, 2).Merge(DTI.Cell(6, 5));
DTI.Cell(6, 3).Merge(DTI.Cell(6, 6));
DTI.Cell(5, 1).Merge(DTI.Cell(6, 1));
DTI.Cell(7, 1).Merge(DTI.Cell(7, 9));
DTI.Cell(8, 1).Merge(DTI.Cell(8, 9));
DTI.Cell(9, 1).Merge(DTI.Cell(9, 3));
DTI.Cell(9, 2).Merge(DTI.Cell(9, 3));
DTI.Cell(9, 3).Merge(DTI.Cell(9, 4));
DTI.Cell(9, 4).Merge(DTI.Cell(9, 5));
DTI.Cell(10, 1).Merge(DTI.Cell(10, 9));
DTI.Cell(11, 5).Merge(DTI.Cell(11, 9));
DTI.Cell(12, 5).Merge(DTI.Cell(12, 9));
DTI.Cell(11, 1).Merge(DTI.Cell(12, 4));

Selection.Start = Content.end;
Selection.TypeParagraph;
Selection.Text = '主管院长签字：            年    月    日';
Paragraphformat.Alignment = 'wdAlignParagraphRight';
Selection.MoveDown;

% 写入表格内容
DTI.Cell(1,1).Range.Text = '课程名称';
DTI.Cell(1,3).Range.Text = '课程号';
DTI.Cell(1,5).Range.Text = '任课教师学院';
DTI.Cell(1,7).Range.Text = '任课教师';
DTI.Cell(2,1).Range.Text = '授课班级';
DTI.Cell(2,3).Range.Text = '考试日期';
DTI.Cell(2,5).Range.Text = '应考人数';
DTI.Cell(2,7).Range.Text = '实考人数';
DTI.Cell(3,1).Range.Text = '出卷方式';
DTI.Cell(3,3).Range.Text = '阅卷方式';
DTI.Cell(3,5).Range.Text = '选用试卷A/B';
DTI.Cell(3,7).Range.Text = '考试时间';
DTI.Cell(4,1).Range.Text = '考试方式';
DTI.Cell(4,3).Range.Text = '平均分';
DTI.Cell(4,5).Range.Text = '不及格人数';
DTI.Cell(4,7).Range.Text = '及格率';
DTI.Cell(5,1).Range.Text = '成绩分布';
DTI.Cell(5,2).Range.Text = '90分以上      人占        %';
DTI.Cell(5,3).Range.Text = '80---89分        人占        %';
DTI.Cell(6,2).Range.Text = '70--79分      人占        %';
DTI.Cell(6,3).Range.Text = '60---69分        人占        %';
DTI.Cell(7,1).Range.Text = ['试卷分析（含是否符合教学大纲、难度、知识覆'...
    '盖面、班级分数分布分析、学生答题存在的共性问题与知识掌握情况、教学中'...
    '存在的问题及改进措施等内容）'];
DTI.Cell(7,1).Range.ParagraphFormat.Alignment = 'wdAlignParagraphLeft';
DTI.Cell(9,2).Range.Text = '签字 :';
DTI.Cell(9,4).Range.Text = '年    月    日';
DTI.Cell(10,1).Range.Text = '教研室审阅意见：';
DTI.Cell(10,1).Range.ParagraphFormat.Alignment = 'wdAlignParagraphLeft';
DTI.Cell(10,1).VerticalAlignment = 'wdCellAlignVerticalTop';
DTI.Cell(11,2).Range.Text = '教研室主任（签字）:          年    月    日';
DTI.Cell(11,2).Range.ParagraphFormat.Alignment = 'wdAlignParagraphLeft';
DTI.Cell(8,1).Range.ParagraphFormat.Alignment = 'wdAlignParagraphLeft';
DTI.Cell(8,1).VerticalAlignment = 'wdCellAlignVerticalTop';
DTI.Cell(9,2).Borders.Item(2).LineStyle = 'wdLineStyleNone';
DTI.Cell(9,2).Borders.Item(4).LineStyle = 'wdLineStyleNone';
DTI.Cell(9,3).Borders.Item(4).LineStyle = 'wdLineStyleNone';
DTI.Cell(11,1).Borders.Item(4).LineStyle = 'wdLineStyleNone';

% 如果当前工作文档中有图形存在，通过循环将图形全部删除
Shape = Document.Shapes;
ShapeCount = Shape.Count;
if ShapeCount ~= 0;
    for i = 1:ShapeCount;
        Shape.Item(1).Delete;
    end;
end;

% 产生正态分布随机数，画直方图，并设置图形属性
zft = figure('units','normalized','position',...
 [0.280469 0.553385 0.428906 0.251302],'visible','off');
set(gca,'position',[0.1 0.2 0.85 0.75]);
rng('default');
data = normrnd(75,6,1000,1);
hist(data);
grid on;
xlabel('考试成绩');
ylabel('人数');

% 将图形复制到粘贴板
hgexport(zft, '-clipboard');
delete(zft);

% 将图形粘贴到当前文档里
% Selection.Range.PasteSpecial;
DTI.Cell(8,1).Range.Paragraphs.Item(1).Range.PasteSpecial;
Shape.Item(1).WrapFormat.Type = 3;
Shape.Item(1).ZOrder('msoBringInFrontOfText');

Document.ActiveWindow.ActivePane.View.Type = 'wdPrintView';
Document.Save;