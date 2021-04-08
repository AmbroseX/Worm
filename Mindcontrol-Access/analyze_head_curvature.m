function analyze_head_curvature(mcd)

if exist('istart', 'var')
        answer = inputdlg({'Start frame', 'End frame'}, 'Cancel to clear previous', 1, ...
            {num2str(istart),num2str(iend)});
    else
        answer = inputdlg({'Start frame', 'End frame'}, '', 1);
end
    
if isempty(answer)
    answer = inputdlg({'Start frame', 'End frame'}, '', 1);
end
    
istart = str2num(answer{1}); %start frame
iend = str2num(answer{2});   %end frame
numframes=iend-istart+1;



T1=mcd(istart).TimeElapsed;
figure;
title('curvature of the head');

for j=2:numframes
    i=istart+(j-1);
    t=mcd(i).TimeElapsed;
    if ~mcd(i).DLPisOn
        plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[mcd(i-1).HeadCurv*100,mcd(i).HeadCurv*100],'k','Linewidth',2);
        hold on;
    else
        if mcd(i).GreenLaser>0
            plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[mcd(i-1).HeadCurv*100,mcd(i).HeadCurv*100],'g','Linewidth',2);
        else
            plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[mcd(i-1).HeadCurv*100,mcd(i).HeadCurv*100],'color',[0.8 1 0.8],'Linewidth',2);
        end
        hold on;
    end
end


figure;
title('Derivative of the curvature near the head');


diff_curv=zeros(numframes,1);
for j=1:numframes
    i=istart+(j-1);
    diff_curv(j)=mcd(i).HeadCurvDeriv;
end

diff_curv_smooth=smooth(diff_curv,10);

for j=2:numframes
    i=istart+(j-1);
    t=mcd(i).TimeElapsed;
    if ~mcd(i).DLPisOn
        %plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[mcd(i-1).HeadCurvDeriv*100,mcd(i).HeadCurvDeriv*100],'k','Linewidth',2);
        plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[diff_curv_smooth(j-1)*100,diff_curv_smooth(j)*100],'k','Linewidth',2);
        hold on;
    else
        if mcd(i).GreenLaser>0
            %plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[mcd(i-1).HeadCurvDeriv*100,mcd(i).HeadCurvDeriv*100],'g','Linewidth',2);
            plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[diff_curv_smooth(j-1)*100,diff_curv_smooth(j)*100],'g','Linewidth',2);
        else
            %plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[mcd(i-1).HeadCurvDeriv*100,mcd(i).HeadCurvDeriv*100],'color',[0.8 1 0.8],'Linewidth',2);
            plot([mcd(i-1).TimeElapsed,mcd(i).TimeElapsed]-T1,[diff_curv_smooth(j-1)*100,diff_curv_smooth(j)*100],'color',[0.8 1 0.8],'Linewidth',2);
        end
        hold on;
    end
end
    








