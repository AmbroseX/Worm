function [se,status,pklcs]=draw_number_peaks
close all
answer=inputdlg({'which sequence do you want to draw','and what status?','the segment you wanna draw'},'',1,...
    {'','2','80'});

se=str2num(answer{1});
status=str2num(answer{2});
segment=str2num(answer{3});
%angle_smt = evalin('base', 't_w_angle_smt');
%peaks=evalin('base', 't_w_peaks');
%t_w_angle=evalin('base', 't_w_angledatazerooffset1');
minimalperiod = evalin('base', 'minimalperiod');
fps = evalin('base', 'fps');
[worm_pks,worm_pklcs] = findpeaks(smooth(smooth(t_w_angle{se,1}{status,1}(:,segment),10,'sgolay'),5));
findpeaks(smooth(smooth(t_w_angle{se,1}{status,1}(:,segment),10,'sgolay'),5));

pklcs=[worm_pks,worm_pklcs];
numpks=length(worm_pklcs);
pks_k=0;
         for pks_j=2:numpks
             
             if pklcs(pks_j-pks_k,2)-pklcs(pks_j-pks_k-1,2)<minimalperiod*fps % get rid of fulse peaks
                 pklcs(pks_j-pks_k,:)=[];
                 pks_k=pks_k+1;
             end
         end
set(gcf,'position',[1200 500 600 400]);
hold on;
text(pklcs(:,2)-5,pklcs(:,1)+0.05,...
    num2str((1:numel(pklcs(:,1)))'));
switch status
    case 1
        illu_status='peaks of head, before DLPon';
    case 2
        illu_status='peaks of head, during DLPon';
    case 3
        illu_status='peaks of tail, during DLPon';
end
title([illu_status,'segment=',answer{3}]);
hold off;
end