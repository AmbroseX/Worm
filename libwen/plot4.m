function h=plot4(x,y,z,c)
    %clf;
    sz = size(x);
    if sz(1)<sz(2)
        x=x';
        y=y';
        z=z';
        c=c';
    end
    if sz(2)==4
        h=surface([x(:,1) x(:,1)],[x(:,2) x(:,2)],[x(:,3) x(:,3)],[x(:,4) x(:,4)],'facecol','no','edgecol','interp','LineWidth',1); colormap jet;
    else
        h=surface([x x],[y y],[z z],[c c],'facecol','no','edgecol','interp','LineWidth',1); colormap jet;
    end
    %hold off;
end