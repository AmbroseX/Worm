
close all;


%disp('If you are analysing part of DLP ON time, just the very number.');
%disp('If you are analysing part of DLP OFF time, then input 100+NUMBER.');
%sequence_extracting=input('which sequence do you want to analyse? \n');
%if sequence_extracting>100
%    sequence_extracting_index=sequence_extracting-100;
%else
%    sequence_extracting_index=sequence_extracting;
%end
%% define the variables
angle_data_selected=cell(cyclenum,1);
answer = inputdlg({'limit for minimal period'}, '', 1);
minimalperiod=str2double(answer{1});
%maxperiod=str2double(answer{2});
t_w=cell(cyclenum,4);
t_w_angle_smt=cell(cyclenum,2);
t_w_angledatazerooffset1=cell(cyclenum,1);
t_w_peaks=cell(cyclenum,1);
%{
h1=figure(1);
set(gcf,'position',[50 500 600 400])
plot(1:100,t_w_amp{sequence_extracting_index,1}(1,:),'k',1:100,t_w_amp{sequence_extracting_index,1}(2,:),'g')
title('body amplitude diagram')
xlabel('body coordinate')
ylabel('amplitude')
legend('before','during');

h2=figure(2);
set(gcf,'position',[620 500 600 400]);
title('click start and end frames to analyze phase velocity and angle of attacks');
%}


%% select frame of start and end point to be analysed
for sequence_extracting_index=1:cyclenum
    % select form the kymograph
    %imagesc(t_w_curvdatafiltered{sequence_extracting_index,1}(:,:)); colormap(cmap); colorbar; caxis([-10 10]); hold on;
    % t1b t2b t1d t2d t1a t2a

    %t1 = ginput(1);
    shortafterflag=0;
    if t_w_illu_time(sequence_extracting_index,1)-300>0
        t1=t_w_illu_time(sequence_extracting_index,1)-300;
    else
        t1=t_w_illu_time(sequence_extracting_index,1)-100;
    end
    if t_w_illu_time(sequence_extracting_index,2)+300<=numframes_total
        t2=t_w_illu_time(sequence_extracting_index,2)+300;
    else
        t2=t_w_illu_time(sequence_extracting_index,2)+50;
        shortafterflag=1;
    end
    t1b=t1-t1+1;
    t2b=t_w_illu_time(sequence_extracting_index,1)-1-t1+1;
    t1d=t_w_illu_time(sequence_extracting_index,1)-t1+1;
    t2d=t_w_illu_time(sequence_extracting_index,2)-t1+1;
    t1a=t_w_illu_time(sequence_extracting_index,2)+1-t1+1;
    t2a=t2-t1+1;
    
    %plot( [1 numcurvpts],[t1b t1b], '-w');
    %t2 = ginput(1);
    %plot( [1 numcurvpts],[t2b t2b], '-w');
    %plot( [1 numcurvpts],[t1d t1d], '-b');
    %plot( [1 numcurvpts],[t1d t1d], '-b');
    %plot( [1 numcurvpts],[t1d t1d], '-g');
    %plot( [1 numcurvpts],[t1d t1d], '-g');

    t_start_auto=t1d;
    t_end_auto=t2d;

    fps=(t_end_auto-t_start_auto)/(time_auto{sequence_extracting_index,1}(t_end_auto,1)-time_auto{sequence_extracting_index,1}(t_start_auto,1));

    c2 = t_w_curvdatafiltered{sequence_extracting_index,1}(t_start_auto:t_end_auto,:) > 0;
    %figure; imagesc(c2(:,:));  hold on;
    %colormap(jet);
    %title('filtered, binary');

    maskhead = 0.1;
    masktail = 0.1;
    minimum_fraction_for_fit = 0.8;

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
    %hold on;

    for n=1:numlab
        c5 = (c4 == n);
        [y, x] = find(c5);

        yshift = 3;
        yshifted = ceil(1+0.5*(1+sign(y-yshift)) .* (y-yshift-1))+t_start_auto;
        curvshift = zeros(size(x));

        for jj=1:length(x)
            curvshift(jj) = t_w_curdata{sequence_extracting_index,1}(yshifted(jj), x(jj));
        end

        tmp = x;
        x=x(logical((tmp>=maskhead * numcurvpts) .* (tmp<=(1-masktail) * numcurvpts)));
        y=y(logical((tmp>=maskhead * numcurvpts) .* (tmp<=(1-masktail) * numcurvpts)));

        if max(x) - min(x) >=  (1-maskhead-masktail)*minimum_fraction_for_fit*numcurvpts

            [p,S] = polyfit(x,y,1);


            if S.normr < normrthresh
                %if mean(curvshift) > 0
                %    plotcol = '-g';
                %else
                %    plotcol = '--g';
                %end
                %plot(polyval(p,[1:numcurvpts]), plotcol); hold on;        
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
                    %text(xpos,ypos,num2str([numcycles2 p(1)]), 'Color', 'white'); hold on;
            end
        end

    end

    %plot( [0.5+maskhead*numcurvpts 0.5+maskhead*numcurvpts],[1 numframes], ':k');
    %plot( [0.5+ (1-masktail)*numcurvpts 0.5+ (1-masktail)*numcurvpts],[1 numframes], ':k');


    %title(strcat(num2str(sum(curvsigndatatmp)),' positive, ', num2str(numcycles2-sum(curvsigndatatmp)),...
    %                                ' negative.  Click on bad fits, press return'));
    %{
    badfits = ginput;

    epsilon = 4;  % how close to fit you need to click
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
    %}
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

    %imagesc(c2); hold on;

    %plot( [0.5+maskhead*numcurvpts 0.5+maskhead*numcurvpts],[1 numframes], ':k');
    %plot( [0.5+ (1-masktail)*numcurvpts 0.5+ (1-masktail)*numcurvpts],[1 numframes], ':k');
   %{
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
    %}
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
    %}
    %% calculate the cycles of wave propagated

    %here wavevelocity represents the the time needs for a complete phase for a
    %whole body length

     angle_data_selected{1,1} = unwrap(unwrap(t_w_angledata{sequence_extracting_index,1}(t1b:t2b,:)),[],2);
     angle_data_selected{2,1} = unwrap(unwrap(t_w_angledata{sequence_extracting_index,1}(t1d:t2d,:)),[],2);
     angle_data_selected{3,1} = unwrap(unwrap(t_w_angledata{sequence_extracting_index,1}(t1a:t2a,:)),[],2);
     angle_data_zerooffset=cell(3,1);
     angle_data_zerooffset1=cell(3,1);
     period_s=zeros(3,1);
     frequency=zeros(3,1);
     wavevelocity=zeros(3,1);
     wavelength=zeros(3,1);
     angle_smth=cell(3);
     angle_smthp20=cell(3,1);
     if ~shortafterflag
         illustatus=3;
     else
         illustatus=2;
     end
     for i=1:illustatus %each i donotes different illumination state
         angle_data_avg1 = mean(angle_data_selected{i,1}(:,:),2);
         angle_data_zerooffset{i,1} = angle_data_selected{i,1} - repmat(angle_data_avg1, [1 size(angle_data_selected{i,1},2)]);
         angle_data_zerooffset{i,1} = unwrap(unwrap(angle_data_zerooffset{i,1},[],2));
         angle_data_zerooffset1{i,1} = angle_data_zerooffset{i,1}-mean(angle_data_zerooffset{i,1}(:));
        % calculat period three different illlumination state by segment 20 of body
         angle_smthp20{i}=smooth(smooth(angle_data_zerooffset1{i,1}(:,20),10,'sgolay'),5);
        %angle_smth data has 3-by-3 dimention, column means illu stats: before, during, after; row means bodysegment
        %15,50,85.
        for j=1:3
            angle_smth{i,j}=smooth(smooth(angle_data_zerooffset1{i,1}(:,15+35*(j-1)),10,'sgolay'),5);
        end
        %{
         figure;clf;
         subplot(121);
         imagesc(t_w_angledatazerooffset1); hold on; colorbar;
         title({'worm angle, zero offset'}, 'Interpreter', 'None');
         xlabel('body coord');
         ylabel('frame #');
         subplot(122);
         imagesc(t_w_angledatazerooffset1>0); hold on; colorbar;
         title({'worm angle, zero offset'}, 'Interpreter', 'None');
         xlabel('body coord');
         ylabel('frame #');
          %}
         [worm_pks,worm_pklcs]=findpeaks(angle_smthp20{i});
         [worm_pksI,worm_pklcsI]=findpeaks(-angle_smthp20{i});
         pklcs=[worm_pks,worm_pklcs];
         pklcsI=[worm_pksI,worm_pklcsI];
         numpks=length(worm_pklcs);
         numpksI=length(worm_pklcsI);
         k=0;
         for j=2:numpks
             
             if pklcs(j-k,2)-pklcs(j-k-1,2)<minimalperiod*fps || pklcs(j-k,1)<0 % get rid of fulse peaks
                 pklcs(j-k,:)=[];
                 k=k+1;
             end
         end
         k=0;
         for j=2:numpksI
             if pklcsI(j-k,2)-pklcsI(j-k-1,2)<minimalperiod*fps || pklcsI(j-k,1)<0 % get rid of fulse peaks
                 pklcsI(j-k,:)=0;
                 k=k+1;
             end
         end
         % stop counitng at reversal dirction
         offflag=0;
         if size(pklcs,1)>1
            for j=2:length(pklcs)
                for k=1:length(pklcsI)
                    if -pklcsI((pklcs(j-1,1)<pklcsI(k,1)<pklcs(j,1)),1)>0
                        offflag=1;
                        break
                    end
                end
                if offflag==1
                    break
                end
            end
         else
             j=1;
         end
         % number of rounds suubstrates the wrong points
         
             rounds_by_frame=(pklcs(j,2)-pklcs(1,2))/(j-1);
         
         period_s(i) = rounds_by_frame/fps;
         frequency(i)=1/period_s(i);
         wavevelocity(i)= mean(1./slopedata)/numcurvpts*fps;
         wavelength(i) = wavevelocity(i) * period_s(i);
         t_w_peak{i,1}=pklcs;
         t_w_peak{i,2}=pklcsI;
         t_w_peak{i,3}=rounds_by_frame;
         t_w_peak{i,4}=j;

     end
     t_w{sequence_extracting_index,1}=period_s;
     t_w{sequence_extracting_index,2}=wavevelocity;
     t_w{sequence_extracting_index,3}=wavelength;
     t_w{sequence_extracting_index,4}=frequency;
     t_w_angle_smt{sequence_extracting_index,1}=angle_smthp20;
     t_w_angle_smt{sequence_extracting_index,2}=angle_smth;
     t_w_angledatazerooffset1{sequence_extracting_index,1}=angle_data_zerooffset1;
     t_w_peaks{sequence_extracting_index,1}=t_w_peak;
     
     clear angle_data_zerooffset curvshift t_w_peak;
     

   
    %if sequence_extracting>100
     %   t_w_meanperiod_DLPoff{sequence_extracting_index,1}=period_s;
      %  t_w_meanperiod_DLPoff{sequence_extracting_index,2}=wavevelocity;
       % t_w_meanperiod_DLPoff{sequence_extracting_index,3}=wavelength;
        %t_w_meanperiod_DLPoff{sequence_extracting_index,4}=sqrt(sin2_theta);

   % end
end
t_w_mat=zeros(size(t_w,1),size(t_w,2)*2);
for j=1:size(t_w,2)
   for i=1:size(t_w,1)
        for k=1:2
            if rem(k,2)
                t_w_mat(i,j)=t_w{i,j}(k,1);
            else
                t_w_mat(i,j+size(t_w,2))=t_w{i,j}(k,1);  
            end
        end
    end
end
%vidfile=strcat(fname(1:end-5),'_HUDS','.avi');

%stagepos=trackShmutz(vidfile,istart,iend);

%TrackWormRealSpace;



close all
