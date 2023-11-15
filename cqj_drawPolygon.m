function [v,pgon] = cqj_drawPolygon()
i = 0;
    grid on,hold on,axis equal
%     axis([-2.4 2.4 -2 2]);
    while 1
        i = i+1;
        try
            [x,y,button] = ginput(1);

        catch
            return;
        end

        v(:,i) = [x;y];
        if i== 1
            pgon = plot(v(1,:),v(2,:),'b-o','MarkerSize',8,'MarkerFaceColor','b');
        else
            set(pgon,'xdata',v(1,:),'ydata',v(2,:));
        end
        if button ~= 1
            set(pgon,'xdata',[v(1,:) v(1)],'YData',[v(2,:) v(2)]);
            break;
        end
    end
end

