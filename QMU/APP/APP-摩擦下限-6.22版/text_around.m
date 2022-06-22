function text_around(x,y,txt,theta,fontsize)
    if nargin==4
        fontsize=10;
    end
    section=mod(theta+pi/12,2*pi);
    if section>pi+pi/6
        %�϶���
        if section>1.5*pi+pi/6
            %�����
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','left','Fontsize',fontsize);
        elseif section>1.5*pi
            %�ж���
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','center','Fontsize',fontsize);
        else
            %�Ҷ���
            text(x,y,txt,'VerticalAlignment','cap','HorizontalAlignment','right','Fontsize',fontsize);
        end
    elseif section>pi
        %�С��Ҷ���
        text(x,y,txt,'VerticalAlignment','middle','HorizontalAlignment','right','Fontsize',fontsize);
    elseif section>pi/6
        %�¶���
        if section>0.5*pi+pi/6
            %�Ҷ���
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',fontsize);
        elseif section>0.5*pi
            %�ж���
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','center','Fontsize',fontsize);
        else
            %�����
            text(x,y,txt,'VerticalAlignment','bottom','HorizontalAlignment','left','Fontsize',fontsize);
        end
    else
        %�С������
        text(x,y,txt,'VerticalAlignment','middle','HorizontalAlignment','left','Fontsize',fontsize);
    end
end

