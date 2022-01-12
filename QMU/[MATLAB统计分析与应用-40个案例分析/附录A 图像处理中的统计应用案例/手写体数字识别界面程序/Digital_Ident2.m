function Digital_Ident2

%手写体数字识别的界面程序
%贝叶斯判别
OldHandle=findobj( 'Type', 'figure', 'Tag', 'shibie' ) ;
if ishandle( OldHandle )
    close( OldHandle ) ;
end
fig=figure ;
init_shibie(fig) ; 
    
%--------------------------------------------------------------------------
function init_shibie(fig)

scnsize=get(0,'screensize');
font_size1=12*sqrt(scnsize(3)*scnsize(4)/(1024*768));
set(fig,'units','normalized','position',[0.15 0.15 0.7 0.6],...
    'menubar','none','name','手写体数字识别',...
    'numbertitle','off','color',[0.925 0.914 0.847],'tag','shibie',...
    'resize','off',...
    'Interruptible','on');
uimenu(gcf,'Label','帮助','Callback',...
    ['nn=sprintf([''            ******手写体数字识别系统******\n\n'' ',...
    ' ''**********************************************************\n\n'' ',...
    ' ''使用前首先点击导入样本图片按钮，\n\n'' ',...
    ' ''导入数字0-9的样本图片，每个字对\n\n'' ',...
    ' ''应5幅图片。选中<选择手绘图片>\n\n'' ',...
    ' ''选项，在坐标系内拖动鼠标写字，\n\n'' ',...
    ' ''然后点击<开始识别>按钮，查看识\n\n'' ',...
    ' ''别结果。也可选中<选择外部>图片\n\n'' ',...
    ' ''选项，选择一幅或多幅外部图片，然\n\n'' ',...
    ' ''后点击<开始识别>按钮，进行识别。\n\n'' ',...
    ' ''<区域选取>按钮可进行识别区域的选\n\n'' ',...
    ' ''择，若未选择区域，则按默认操作。'']);'...
    'h=msgbox(nn,''帮助'',''replace'');'])
uimenu(gcf,'Label','关于','Callback',...
    ['nn=sprintf([''手写体数字识别系统\n\n'' ',...
    ' ''作者：谢中华\n\n'' ',...
    ' ''QQ：56330069\n\n'' ',...
    ' ''Email：xiezhh@tust.edu.cn\n\n'' ',...
    ' ''制作时间：2009年4月\n\n'' ',...
    ' ''天津科技大学理学院数学系'']);'...
    'h=msgbox(nn,''关于'',''replace'');'])
axes('pos',[0.06 0.1 0.6 0.85],'tag','axes1',...
    'Box','on','ButtonDownFcn',{@mouse_xiezh 'start'});
axis image
backcolor=get(gcf,'Color');
uicontrol(gcf,'style','text','units','normalized','pos',[0.05 0.01 0.5 0.05],...
  'string','','fontsize',font_size1,'fontweight','bold',...
  'fontunits','normalized','Hor','left',...
  'ForegroundColor',[0 0 0],'bac',backcolor,...
  'tag','pictureselectioninfo',...
  'visible','off');
panel1=uipanel(fig,'units','normalized','pos',[0.704 0.01 0.286 0.95]);
panel2=uipanel(fig,'units','normalized','pos',[0.712 0.79 0.27 0.15],...
    'title','输出识别结果',...
    'fontsize',font_size1,'fontweight','bold');
uicontrol(panel2,'style','edit','units','normalized',...
     'pos',[0.1 0.15 0.8 0.7],...
     'fontsize',font_size1,'fontunits','normalized',...
     'fontweight','bold','string','',...
     'tag','edit1','backgroundcolor',[1 1 1]);
 uicontrol(panel1,'style','push','units','normalized',...
     'pos',[0.1 0.71 0.8 0.08],...
     'string','导入样本图片',...
     'fontsize',font_size1,'fontweight','bold',...
     'fontunits','normalized',...
     'tag','ImportPicture',...
     'callback',@SamplesPicture)
 
 panel3=uipanel(fig,'units','normalized','pos',[0.712 0.45 0.27 0.2],...
    'title','选择要识别的图片',...
    'fontsize',font_size1,'fontweight','bold');
 uicontrol(panel3,'style','radiobutton','units','normalized',...
    'pos',[0.12 0.1 0.75 0.35],'tag','ExternalImage',...
    'value',0,'string','选择外部图片',...
    'fontsize',font_size1,'fontweight','bold',...
    'fontunits','normalized','callback',@pictureselection_callback2);
 uicontrol(panel3,'style','radiobutton','units','normalized',...
    'pos',[0.12 0.55 0.75 0.35],'tag','HandPaintedPicture',...
    'value',1,'string','选择手绘图片',...
    'fontsize',font_size1,'fontweight','bold',...
    'fontunits','normalized','callback',@pictureselection_callback1);
panel4=uipanel(fig,'units','normalized','pos',[0.712 0.24 0.27 0.18],...
    'title','选择绘图颜色',...
    'fontsize',font_size1,'fontweight','bold');
uicontrol(panel4,'style','popupmenu',...
     'units','normalized',...
     'pos',[0.1 0.3 0.8 0.5],...     
     'string','黑色|蓝色|红色|黄色|绿色|品红|青色',...
     'fontsize',font_size1,'fontunits','normalized',...
     'backgroundcolor',[1 1 1],...
     'tag','huitu_color','value',1);
uicontrol(panel1,'style','push','units','normalized',...
     'pos',[0.1 0.13 0.38 0.08],...
     'string','区域选取',...
     'fontsize',font_size1,'fontweight','bold',...
     'fontunits','normalized',...
     'tag','RegionSelection',...
     'callback',@regionselection_callback) 
uicontrol(panel1,'style','push','units','normalized',...
     'pos',[0.52 0.13 0.38 0.08],...
     'string','开始识别',...
     'fontsize',font_size1,'fontweight','bold',...
     'fontunits','normalized',...
     'tag','BeginIdent',...
     'callback',@Digital_Identification_call)
uicontrol(panel1,'style','push','units','normalized',...
     'pos',[0.1 0.02 0.38 0.08],...
     'string','清除',...
     'fontsize',font_size1,'fontweight','bold',...
     'fontunits','normalized',...
     'tag','RemovePicture',...
     'callback',@clearaxes_callback)
uicontrol(panel1,'style','push','units','normalized',...
     'pos',[0.52 0.02 0.38 0.08],...
     'string','退出',...
     'fontsize',font_size1,'fontweight','bold',...
     'fontunits','normalized',...
     'tag','QuitGui',...
     'callback',['close 手写体数字识别;'])

setappdata(gcf,'xy_point',[]);
handles=guihandles(gcf);
guidata(gcf,handles);

%%-------------------------------------------------------------------------
%鼠标手绘
%%-------------------------------------------------------------------------
function mouse_xiezh(hObject, eventdata, action)

handles=guidata(gcf);
data.color=get(handles.huitu_color,'value');

switch data.color
    case 1
        str1='k';
    case 2
        str1='b';
    case 3
        str1='r';
    case 4
        str1='y';
    case 5
        str1='g';
    case 6
        str1='m';
    case 7
        str1='c';        
end

switch action
    case 'start'        
        point=get(gca,'CurrentPoint'); 
        xy_point=getappdata(gcf,'xy_point');
        xy_point=[xy_point;point(1,1:2)];
        setappdata(gcf,'xy_point',xy_point);
        line(point(1,1),point(1,2),'clipping','on',...
            'erasemode','background','linewidth',3,'color',str1);
        data_point=point(1,1:2);        
        setappdata(gcf,'data_point',data_point);
        set(gcbf,'WindowButtonMotionFcn',{@mouse_xiezh 'move'});
        set(gcbf,'WindowButtonUpFcn',{@mouse_xiezh 'stop'});        
    case 'move'
        data_point=getappdata(gcf,'data_point');
        point=get(gca,'CurrentPoint'); 
        xy_point=getappdata(gcf,'xy_point');
        xy_point=[xy_point;point(1,1:2)];
        setappdata(gcf,'xy_point',xy_point);
        line([data_point(1),point(1,1)],[data_point(2),point(1,2)],...
            'clipping','on','erasemode','background',...
            'linewidth',3,'color',str1); 
        data_point=point(1,1:2);
        setappdata(gcf,'data_point',data_point);
    case 'stop'
        set(gcbf,'WindowButtonMotionFcn','');
        set(gcbf,'WindowButtonUpFcn','');
end

%%-------------------------------------------------------------------------
%导入样本图片
%%-------------------------------------------------------------------------
function SamplesPicture(hObject, eventdata)

geshi={'*.jpg','JPEG image (*.jpg)';...
       '*.bmp','Bitmap image (*.bmp)';...
       '*.*','All Files (*.*)'};
[FileName FilePath]=uigetfile(geshi,'导入外部数据','*.jpg','MultiSelect','on');
if ~isequal([FileName,FilePath],[0,0]);
    FileFullName=strcat(FilePath,FileName);
%     if  ~ischar(FileFullName)
%         FileFullName=FileFullName([2:end 1])';
%     end    
else
    return;
end

n = length(FileFullName);
I = zeros(50);
BW = zeros(16);
training = zeros(1,256);
group = [];
for i=1:n
    I=imread(FileFullName{i});
    %if ndims(I)==3
    %    I=255-I(:,:,1);
    %end
    I=im2bw(255-I,0.4);
    [y,x]=find(I==1);
    BW=I(min(y):max(y),min(x):max(x));
    BW=imresize(BW,[16 16]);
    training(i,:) = BW(:)';    
    [pathstr,namestr,ext]=fileparts(FileName{i});
    group = [group;str2num(namestr(4))];
end
ObjBayes = NaiveBayes.fit(training,group,'Distribution','mn');
save shuzi_ObjBayes.mat ObjBayes
msgbox('样本图片导入完毕','消息','replace');


%%-------------------------------------------------------------------------
%数字识别
%%-------------------------------------------------------------------------
function Digital_Identification_call(hObject, eventdata)

if exist('shuzi_ObjBayes.mat')
    load shuzi_ObjBayes.mat
    handles=guidata(gcf);
    HandPainted_radio=get(handles.HandPaintedPicture,'value');
    External_radio=get(handles.ExternalImage,'value'); 
    I=[];
    if HandPainted_radio==1
        xy_point=getappdata(gcf,'xy_point');
        if ~isempty(xy_point)
            x_point=xy_point(:,1);
            y_point=xy_point(:,2);
            tag_select=getappdata(gcf,'tag_select');
            if tag_select==1
                xlim1=xlim;
                ylim1=ylim;
                setappdata(gcf,'tag_select',0);
                if min(x_point) > xlim1(2) | max(x_point) < xlim1(1) | min(y_point) > ylim1(2) | max(y_point) < ylim1(1) 
                    warndlg('您选择的区域里可能还没有写数字吧！！！','警告','replace');                    
                    return;                    
                else
                    new_fig=figure('color',[1 1 1],'visible','off','position',[0.5 0.5 0.2 0.2]);
                    new_axes=copyobj(handles.axes1,new_fig);
                    set(new_axes,'position',[0.1 0.1 0.8 0.8]);                     
                end
            else                
                new_fig=figure('color',[1 1 1],'visible','off','position',[0.5 0.5 0.2 0.2]);
                new_axes=copyobj(handles.axes1,new_fig);
                set(new_axes,'position',[0.1 0.1 0.8 0.8]);            
                xrange=range(x_point);
                yrange=range(y_point);
                if xrange>=yrange
                    h=xrange-yrange;
                    axis([min(x_point) max(x_point) min(y_point)-h/2 max(y_point)+h/2]);
                else
                    h=yrange-xrange;
                    axis([min(x_point)-h/2 max(x_point)+h/2 min(y_point) max(y_point)]);
                end
            end    
            axis off;
            set(gcf,'PaperPosition',[0 0 0.85 0.85],'PaperUnits','normalized')
            saveas(gca,[pwd '\zi.jpg']);
            delete(new_fig);            
            I=imread([pwd '\zi.jpg']); 
            BW=im2bw(255-I,0.4);
            Digital_Identification(BW,ObjBayes);
            %BW=bwareaopen(BW,10);  %去除干扰点
            %se=strel('disk',15);
            %BW=imdilate(BW,se); %图像膨胀
            setappdata(gcf,'xy_point',[]);
        else
            warndlg('您可能还没有写数字吧！！！','警告','replace');
        end
    end
    if External_radio==1
        FileFullName=getappdata(gcf,'FileFullName');
        FileName=getappdata(gcf,'FileName');
        if ~isempty(FileFullName)
            picturenumber=size(FileFullName,1);
            setappdata(gcf,'FileFullName',[]);
            if picturenumber==1;
                BW=getappdata(gcf,'BW');
                %se=strel('disk',3);
                %BW=imdilate(BW,se); %图像膨胀
                Digital_Identification(BW,ObjBayes);
            else
                h=waitbar(0,'正在识别，请稍候...');
                for i=1:picturenumber
                    I=imread(FileFullName{i});
                    BW=im2bw(255-I,0.4);
                    %se=strel('disk',3);
                    %BW=imdilate(BW,se); %图像膨胀
                    result{i}=Digital_Identification(BW,ObjBayes);
                    waitbar(i/picturenumber,h,['当前已完成',num2str(100*i/picturenumber),'%...']);
                    drawnow;
                end                
                result=result';
                delete(h);
                xlswrite('识别结果.xls',{'文件名','识别结果'},'A1:B1');
                xlswrite('识别结果.xls',FileName,['A2:A' num2str(picturenumber+1)]);
                xlswrite('识别结果.xls',result,['B2:B' num2str(picturenumber+1)]);
                set(handles.edit1,'string','识别结果已导入Excel文件');
            end            
        else
            warndlg('您可能还没有选择图片！！！','警告','replace');
        end
    end
else
    warndlg('您可能还没有导入样本图片！！！','警告','replace');
end

%%-------------------------------------------------------------------------
%清空当前坐标系
%%-------------------------------------------------------------------------
function clearaxes_callback(hObject, eventdata)

handles=guidata(gcf);
HandPainted_radio=get(handles.HandPaintedPicture,'value');
set(handles.edit1,'string','');
setappdata(gcf,'xy_point',[]);
cla;
axis([0 1 0 1]);
if HandPainted_radio==1
    set(handles.axes1,'visible','on',...
      'TickDir','in','YDir','normal','layer','bottom', ...
      'ButtonDownFcn',{@mouse_xiezh 'start'});
else
    set(handles.axes1,'visible','on',...
      'TickDir','in','YDir','normal','layer','bottom', ...
      'ButtonDownFcn',[]);
end
set(handles.pictureselectioninfo,'visible','off');

%%-------------------------------------------------------------------------
%选取识别区域
%%-------------------------------------------------------------------------
function regionselection_callback(hObject, eventdata)

handles=guidata(gcf);
HandlesGroup=[handles.ImportPicture;...
    handles.ExternalImage;...
    handles.HandPaintedPicture;...
    handles.RegionSelection;...
    handles.BeginIdent;...
    handles.RemovePicture;...
    handles.QuitGui];
HandPainted_radio=get(handles.HandPaintedPicture,'value');
External_radio=get(handles.ExternalImage,'value'); 
set(HandlesGroup,'enable','off');
if HandPainted_radio==1
    xy_point=getappdata(gcf,'xy_point');
    if ~isempty(xy_point)
        rectsel=getrect(handles.axes1);
        if ~any(rectsel(3:4)==0)
            axis([rectsel(1),rectsel(1)+rectsel(3),rectsel(2),rectsel(2)+rectsel(4)]);
            setappdata(gcf,'tag_select',1);
        end
    end
end
if External_radio==1
    FileFullName=getappdata(gcf,'FileFullName');
    if ~isempty(FileFullName)
        picturenumber=size(FileFullName,1); 
        if picturenumber==1;
            BW=getappdata(gcf,'BW');
            BW=imcrop(BW);
            setappdata(gcf,'BW',BW);
            imshow(BW);
        else
            warndlg('您选择了多幅图片，不能再选择识别区域！！！','警告','replace');
        end
    end
end
set(HandlesGroup,'enable','on');

%%-------------------------------------------------------------------------
%选择手绘图片
%%-------------------------------------------------------------------------
function pictureselection_callback1(hObject, eventdata)

handles=guidata(gcf);
set(handles.axes1,'visible','on',...
  'TickDir','in','YDir','normal','layer','bottom', ...
  'ButtonDownFcn',{@mouse_xiezh 'start'});
set(handles.ExternalImage,'value',0);
set(handles.HandPaintedPicture,'value',1);

%%-------------------------------------------------------------------------
%选择外部图片
%%-------------------------------------------------------------------------
function pictureselection_callback2(hObject, eventdata)

handles=guidata(gcf);
set(handles.axes1,'ButtonDownFcn',[]);
set(handles.ExternalImage,'value',1);
set(handles.HandPaintedPicture,'value',0);
FileFullName=[];
geshi={'*.jpg','JPEG image (*.jpg)';...
   '*.bmp','Bitmap image (*.bmp)';...
   '*.*','All Files (*.*)'};
[FileName FilePath]=uigetfile(geshi,'导入外部数据','*.jpg','MultiSelect','on');
if ~isequal([FileName,FilePath],[0,0]);
    FileFullName=strcat(FilePath,FileName);
    if  ischar(FileFullName)
        I=imread(FileFullName);
    else
        FileFullName=FileFullName([2:end 1])';
        FileName=FileName([2:end 1])';
        I=imread(FileFullName{1});
    end    
    BW=im2bw(255-I,0.4);
    setappdata(gcf,'BW',BW);
    %BW=bwareaopen(BW,10);  %去除干扰点
    %se=strel('disk',3);
    %BW=imdilate(BW,se); %图像膨胀
    imshow(BW);
    picturenumber=size(FileFullName,1); 
    if picturenumber>1
        set(handles.pictureselectioninfo,...
          'string',['您选择了 ' num2str(picturenumber) ' 幅图片，这里只显示第一幅。'],...
          'visible','on');
    end
end
setappdata(gcf,'FileFullName',FileFullName);
setappdata(gcf,'FileName',FileName);

%%-------------------------------------------------------------------------
%公用识别函数
%%-------------------------------------------------------------------------
function result_shu = Digital_Identification(BW,ObjBayes)

handles=guidata(gcf);
result_shu=[]; 
if any(any(BW))
    [y,x]=find(BW==1);
    BW=BW(min(y):max(y),min(x):max(x));
    BW=imresize(BW,[16 16]);
    %imwrite(BW,'zi.bmp');
    %figure
    if nargout==0
        imshow(BW)
    end
    BW=double(BW(:)');
    result_shu=ObjBayes.predict(BW);
%     if tag<=0.07 & range(juli)<=3
%         if nargout==0
%             warndlg('您写的是嘛！，拜托您写一个大一点规范点的数字（0-9），不要写汉字和英文来考验我！','警告','replace');
%             return;
%         else
%             result_shu=[]; 
%         end
%     end    
else
    if nargout==0
        warndlg('您识别的可能是一幅没写数字的图片！！！','警告','replace');
    end
end
if nargout==0
    set(handles.edit1,'string',num2str(result_shu));
end