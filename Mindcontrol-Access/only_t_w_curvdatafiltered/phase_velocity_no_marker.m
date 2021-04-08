close all;
%if ~exist('t_w_wave_v','var')
%    t_w_meanperiod_DLPon=zeros(cyclenum,2);
%end
%if ~exist('t_w_meanperiod_DLPoff','var')
%    t_w_meanperiod_DLPoff=zeros(cyclenum,2);
%end
%if ~exist('t_w_curvdatafiltered','var')
%    t_w_curvdatafiltered = this_worm_curvdatafiltered;
    %t_w_illu_time = this_worm_illumination_time;
    %t_w_curdata = this_worm_curdata;
%end
answer = inputdlg('worm name','',1);
t_w_curvdatafiltered = evalin('base',[answer{1},'_curvdatafiltered']);
for sequence_extracting = 1:length(t_w_curvdatafiltered)
%if timing_selection == 1
t1 = 1;
t2 = size(t_w_curvdatafiltered{sequence_extracting},1);

%end

t_start = t1;
t_end = t2;

fps=numframes_total/(mcd(numframes_total).TimeElapsed-mcd(1).TimeElapsed);

c2 = t_w_curvdatafiltered{sequence_extracting,1}(t_start:t_end,:) > 0;
figure; imagesc(c2(:,:));  hold on;
colormap(jet);
title('filtered, binary');

maskhead = 0.1;
masktail = 0.2;
minimum_fraction_for_fit = 0.4;

c3 = edge(single(c2),'sobel', 0);

c3(:,round((1-masktail)*numcurvpts):end) = 0;
c3(:,1:round(maskhead*numcurvpts)) = 0;
        
[c4,numlab] = bwlabel(c3);

numcycles2 = 0;

clear slopedata timedata slopedatatmp timedatatmp okdatatmp curvsigndatatmp curvsigndata;

okdatatmp = zeros(numlab, 1);

normrthresh = 220;

%   draw fit limits

%     
hold on;
fitted_line_flags = 0;
for n=1:numlab
    c5 = (c4 == n);
    [y, x] = find(c5);

    yshift = 3;
    yshifted = ceil(1+0.5*(1+sign(y-yshift)) .* (y-yshift-1));
    curvshift = zeros(size(x));
            
    for jj=1:length(x)
        curvshift(jj) = t_w_curvdatafiltered{sequence_extracting,1}(yshifted(jj), x(jj));
    end
    %         
    %         disp(n);
    %         disp('mean(c2shift)');
    %         disp(mean(curvshift));
    %         
            % exclude points near head and tail for fitting
    tmp = x;
    x=x(logical((tmp>=maskhead * numcurvpts) .* (tmp<=(1-masktail) * numcurvpts)));
    y=y(logical((tmp>=maskhead * numcurvpts) .* (tmp<=(1-masktail) * numcurvpts)));

    if max(x) - min(x) >=  (1-maskhead-masktail)*minimum_fraction_for_fit*numcurvpts
        
        [p,S] = polyfit(x,y,1);
    %             plot(polyval(p,[1:numcurvpts]), 'r'); hold on;    
        %         disp([n , S.normr , S.normr < normrthresh]);
    %             S.normr

        if S.normr < normrthresh
            fitted_line_flags = fitted_line_flags+1;
            if mean(curvshift) > 0
                plotcol = '-g';
            else
                plotcol = '--g';
            end
            plot(polyval(p,[1:numcurvpts]), plotcol); hold on;        
            numcycles2 = numcycles2 + 1;
            slopedatatmp(n) = p(1);
            timedatatmp(n) = p(2);
            okdatatmp(n) = 1;
            negshift = (mean(curvshift) > 0);
            curvsigndatatmp(n) = negshift;
            xpos = 5;
            ypos = p(2)-1;
            if p(2)<1
                xpos = numcurvpts/4;
                ypos = 5;
            end
                text(xpos,ypos,num2str([numcycles2 p(1)]), 'Color', 'white'); hold on;
        end
    
    end
   
end
if fitted_line_flags > 0
    plot( [0.5+maskhead*numcurvpts 0.5+maskhead*numcurvpts],[1 numframes], ':k');
    plot( [0.5+ (1-masktail)*numcurvpts 0.5+ (1-masktail)*numcurvpts],[1 numframes], ':k');


    title(strcat(num2str(sum(curvsigndatatmp)),' positive, ', num2str(numcycles2-sum(curvsigndatatmp)),...
                                    ' negative.  Click on bad fits, press return'));

    badfits = ginput;

    epsilon = 8;  % how close to fit you need to click
    for j=1:size(badfits,1)
        for n=1:numlab
            if okdatatmp(n) 
                if abs(timedatatmp(n) + slopedatatmp(n)*badfits(j,1) - badfits(j,2))<epsilon
                    %                 disp(strcat('matches #', num2str(n)));
                    okdatatmp(n) = 0;
                end
            end
        end
    end

    numcycles2 = 0;
    c4b = c4;
    for n=1:numlab
        if okdatatmp(n)
            numcycles2 = numcycles2+1;
            slopedata(numcycles2) = slopedatatmp(n);
            timedata(numcycles2) = timedatatmp(n);
            curvsigndata(numcycles2)=curvsigndatatmp(n);
            c4b(c4b==n) = numcycles2;
        else
            c4b(c4b==n) = 0;
        end
    end

    imagesc(c2); hold on;

    plot( [0.5+maskhead*numcurvpts 0.5+maskhead*numcurvpts],[1 numframes], ':k');
    plot( [0.5+ (1-masktail)*numcurvpts 0.5+ (1-masktail)*numcurvpts],[1 numframes], ':k');

    for n=1:numcycles2
        if curvsigndata(n)
            plotcol = '-g';
        else
            plotcol = '--g';
        end

            plot(polyval([slopedata(n) timedata(n)],[1:numcurvpts]), plotcol); hold on;    
            %             text(5,p(2)-1,num2str([n numcycles2 S.normr]), 'Color', 'white'); hold on;
            %             if p(2)<1 
            %                 text(5,2,num2str([n numcycles2 S.normr]), 'Color', 'white'); hold on;
            %             end
    end

    %if mean(curvsigndata) ~= 0.5
    %    msgbox('Warning: unequal number of positive and negative fits','','error') 
    %end

    %title('Press return to continue');
    %pause;
    %{
    h3=figure(3);
    set(gcf,'position',[1240 500 600 400]);
    imagesc(c2(:,:));

    title('click two points separated in time_auto by N cycles');
    t1 = ginput(1);

    plot( [1 numcurvpts],[t1(2) t1(2)], '-w');
    t2 = ginput(1);

    plot( [1 numcurvpts],[t2(2) t2(2)], '-w');
        %   draw fit limits

    plot( [0.5+maskhead*numcurvpts 0.5+maskhead*numcurvpts],[1 numframes], ':k');
    plot( [0.5+ (1-masktail)*numcurvpts 0.5+ (1-masktail)*numcurvpts],[1 numframes], ':k');

    answer = inputdlg('Enter number of cycles');
    period = abs(t1(2) - t2(2)) / str2num(answer{1});

    title(strcat(answer{1},' cycles, ', num2str(numcycles2), ' fits'), 'Interpreter', 'None');


    for n=1:numcycles2
        if curvsigndata(n)
            plotcol = '-g';
        else
            plotcol = '--g';
        end

        plot(polyval([slopedata(n) timedata(n)],[1:numcurvpts]), plotcol); hold on;    

            %             text(5,p(2)-1,num2str([n numcycles2 S.normr]), 'Color', 'white'); hold on;
            %             if p(2)<1 
            %                 text(5,2,num2str([n numcycles2 S.normr]), 'Color', 'white'); hold on;
            %             end
    end


    period_s = period/fps;
    %fprintf('the %d data are',segseq)
    disp('mean period (sec)')
    disp(period_s);
    disp('frequency(Hz)')
    disp(1/period_s);
    %}
    if exist('slopedata','var')
        wavevelocity = mean(1./slopedata)/numcurvpts*fps;
        clear slopedata timedata curvsigndata
    else
        wavevelocity = 0;
    end
else
    wavevelocity = 0;
end
fprintf('No. %d mean wave velocity (L/s)\n',sequence_extracting)
disp(wavevelocity);
%{
wavelength = wavevelocity * period_s;
disp('wavelength = velocity*period (L)');
disp(wavelength);



 angle_data_selected = unwrap(unwrap(angle_data(t_start:t_end,:),1),2);
 angle_data_avg = imfilter(angle_data_selected, ones(round(period),1)/round(period), 'replicate', 'same');
 angle_data_avg1 = mean(angle_data_avg, 2);
 % figure(104);clf;
 % imagesc(wormangle_avg)
 % figure(105);clf;
 % plot(wormangle_avg(:,1));

 angle_data_zerooffset = angle_data_selected - repmat(angle_data_avg1, [1 size(angle_data_selected,2)]);
 angle_data_zerooffset = unwrap(unwrap(angle_data_zerooffset,2),1);
 figure;clf;
 subplot(121);
 imagesc(angle_data_zerooffset); hold on; colorbar;
 title({'worm angle, zero offset'}, 'Interpreter', 'None');
 xlabel('body coord');
 ylabel('frame #');
 subplot(122);
 imagesc(angle_data_zerooffset>0); hold on; colorbar;
 title({'worm angle, zero offset'}, 'Interpreter', 'None');
 xlabel('body coord');
 ylabel('frame #');
  
 angle_data_zerooffset1 = angle_data_zerooffset- mean(angle_data_zerooffset(:));
 pause(1);
         
         
%        figure(8);subplot(222);
%         plot(angle_data_zerooffset1(:,1:10:end));
%         xlim([0 size(angle_data_zerooffset1,1)]);
%         title({title1; 'wormangle, 1 cycle'}, 'Interpreter', 'None');
%         xlabel('frame #');
%         ylabel('rad');

angle_data_zerooffset1_max = max(angle_data_zerooffset1,[],1);
angle_data_zerooffset1_min = min(angle_data_zerooffset1,[],1);
angle_data_zerooffset1_ampl =angle_data_zerooffset1_max - angle_data_zerooffset1_min;
angle_data_zerooffset1_ampl_by_rms =2* sqrt(2) * std(angle_data_zerooffset1,1);
%         angle_data_zerooffset1_ampl_by_rms =2* sqrt(2) * std(angle_data_zerooffset,1);
%         sin2_theta = mean(sin(angle_data_zerooffset1_ampl/2).^2);
sin2_theta = mean(sin(angle_data_zerooffset1_ampl_by_rms/2).^2);


figure;
title('angle');
plot(angle_data_zerooffset1_max, '-b');hold on; 
plot(-angle_data_zerooffset1_min, '-g');
plot(angle_data_zerooffset1_ampl, '-r');
plot(angle_data_zerooffset1_ampl_by_rms, '-c');
title({ 'wormangle amplitude'});
xlabel('body coord', 'Interpreter', 'None');
xlim([0 100]);
ylabel('rad');

disp('sin_theta');
disp(sqrt(sin2_theta));
disp('mean angle of attack');
disp(mean(angle_data_zerooffset1_ampl_by_rms/2));
disp('standard deviation of head curvature');
disp(std(head_curv(t_start:t_end)));
%}

    %t_w_meanperiod_DLPon(sequence_extracting,1)=period_s;
    %t_w_meanperiod_DLPon(sequence_extracting,2)=1/period_s;
    t_w_wave_v(sequence_extracting,2)=abs(wavevelocity);
    evalc([wormname,'_wave_v = t_w_wave_v']);
%    t_w_meanperiod_DLPon(sequence_extracting,3)=wavelength;
%    t_w_meanperiod_DLPon(sequence_extracting,4)=sqrt(sin2_theta);
%    t_w_meanperiod_DLPon(sequence_extracting,5)=mean(angle_data_zerooffset1_ampl_by_rms/2);
%    t_w_meanperiod_DLPon(sequence_extracting,5)=std(head_curv(t_start:t_end));
end
save([wormname,'_wave_v'],[wormname,'_wave_v']);



