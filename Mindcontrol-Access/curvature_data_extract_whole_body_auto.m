%This version contains the completed functions of analysing worm behavior
%parameter.
%analyse + - 200 frames marked by DLP on signal.(161024, for WEN0056)

disp('Is the worm illuminated?');
illu_status=input('1 for yes, 0 for no\n');
cyclenum=0;
for cycle=1:size(t_w_illu_time)
    if t_w_illu_time(cycle,1)
        cyclenum=cyclenum+1;
    end
end
illu_times=0;
for k=1:size(t_w_illu_time)
    if t_w_illu_time(k,1)
        illu_times=illu_times+1;
    end
end
time_auto=cell(illu_times,1);
frames_beforeDLPon=310
for cycle=1:cyclenum 
    close all
    clear amp curdata curvdatafilterd angle_data
    if t_w_illu_time(cycle,1)<=frames_beforeDLPon
        istart = t_w_illu_time(cycle,1);
    else
        istart = t_w_illu_time(cycle,1)-frames_beforeDLPon;
    end
    if t_w_illu_time(cycle,2)+300>=length(mcd)
        iend = t_w_illu_time(cycle,2)+50;
    else
        iend = t_w_illu_time(cycle,2)+50;
    end
    spline_p = 0.0005;
    flip = 0;

    numframes=iend-istart+1;

    numcurvpts=100;

    proximity = 50;

    curvdata=zeros(numframes,numcurvpts);
    angle_data = zeros(numframes,numcurvpts+1);
    Head_position=mcd(istart).Head;
    Tail_position=mcd(istart).Tail;

    worm_length=0;  %body length in terms of pixels

    mcd2=mcd;

    t1=0;

    j1=0; j2=0;

    Centerline=zeros(numframes,100,2);
    time_auto{cycle,1}=zeros(numframes,1);
    %Display the target regoin in the screen

    for j=1:numframes
        i = istart + (j - 1);
        if ~mcd(i-1).DLPisOn && mcd(i).DLPisOn
            origin_y=mcd(i).IllumRectOrigin(2);
            radius_y=mcd(i).IllumRectRadius(2);
        if origin_y-radius_y<0
            initiation=0;
        else
            initiation=origin_y-radius_y;
        end
            fprintf('\n ************************** \n the target region is [%d,%d].\n **************************\n',initiation,origin_y+radius_y);
            DLPon_frame=i;
        elseif mcd(i-1).DLPisOn && ~mcd(i).DLPisOn
            origin_y=mcd(i).IllumRectOrigin(2);
            radius_y=mcd(i).IllumRectRadius(2);
        if origin_y-radius_y<0
            initiation=0;
        else
            initiation=origin_y-radius_y;
        end
            fprintf(' ************************** \n the target region is [%d,%d].\n **************************\n',initiation,origin_y+radius_y);
            DLPoff_frame=i;
            break;
        end
    end

    for j=1:numframes

        i = istart + (j - 1);

        if (norm(mcd(i).Head-Head_position)> norm(mcd(i).Tail-Head_position)) %%head and tail flips
            if norm(mcd(i).Head-Tail_position)<=proximity && norm(mcd(i).Tail-Head_position)<=proximity  %%if the tip points are identified
                flip=1;
                Head_position=mcd(i).Tail;
                Tail_position=mcd(i).Head;
                %mcd2(i).Head=Head_position;
                %mcd2(i).Tail=Tail_position;
            end
        else
            flip =0;
            Head_position = mcd(i).Head;
            Tail_position = mcd(i).Tail;
        end


        %if  norm(mcd2(i).Head-Head_position)<=proximity && norm(mcd2(i).Tail-Tail_position)<=proximity
        if norm(mcd(i).Head-mcd(i).Tail)>proximity
             centerline=reshape(mcd(i).SegmentedCenterline,2,[]);
             %Head_position=mcd2(i).Head;
             %Tail_position=mcd2(i).Tail;
            if flip
                centerline(1,:)=centerline(1,end:-1:1);
                centerline(2,:)=centerline(2,end:-1:1);
            end
        end


        Centerline(j,:,1)=centerline(1,:);
        Centerline(j,:,2)=centerline(2,:);




        time_auto{cycle,1}(j,1)=mcd(i).TimeElapsed;

        if mcd(i).DLPisOn && ~mcd(i-1).DLPisOn
            t1=time_auto{cycle,1}(j,1);%DLPon timepoint
            w1=t1;
            j1=j;%DLPon frame
            %origin=100-mcd(i).IllumRectOrigin(2);
            %radius=mcd(i).IllumRectRadius(2);
        end

        if ~mcd(i).DLPisOn && mcd(i-1).DLPisOn
            t2=time_auto{cycle,1}(j,1);%DLPoff timepoint
            w2=t2;
            j2=j;%DLPoff frame
        end



        figure (1);
        plot(centerline(1,:),centerline(2,:),'k-');
        hold on; plot(Head_position(1),Head_position(2),'ro');
        hold on; plot(Tail_position(1),Tail_position(2),'bo');

        axis off; axis equal; hold on;
        df = diff(centerline,1,2); 
        t = cumsum([0, sqrt([1 1]*(df.*df))]); 
        worm_length=worm_length+t(end);
        cv = csaps(t,centerline,spline_p);%smooth the curvature

        figure(1);
        fnplt(cv, '-g'); hold off;   

        cv2 =  fnval(cv, t)';
        df2 = diff(cv2,1,1); df2p = df2';

        splen = cumsum([0, sqrt([1 1]*(df2p.*df2p))]);
        cv2i = interp1(splen+.00001*[0:length(splen)-1],cv2, [0:(splen(end)-1)/(numcurvpts+1):(splen(end)-1)]);

        df2 = diff(cv2i,1,1);
        atdf2 =  unwrap(atan2(-df2(:,2), df2(:,1)));
        angle_data(j,:) = atdf2';

        curv = unwrap(diff(atdf2,1)); 
        curvdata(j,:) = curv';



    end
    vedio_sequence=cycle;
    vedio_sequence_s=num2str(vedio_sequence);
    vedio_sequence_s=strcat(wormname,'_',vedio_sequence_s);

    cmap=redgreencmap;
    cmap(:,3)=cmap(:,2);
    cmap(:,2)=0;
    origin=10;
    radius=8;

    worm_length=worm_length/numframes;

    %answer = inputdlg({'time filter', 'body coord filter', 'mean=0, median=1'}, '', 1, {num2str(5), num2str(10), '0'});
    timefilter =5;
    bodyfilter =10;



    h = fspecial('average', [timefilter bodyfilter]);
    curvdatafiltered = imfilter(curvdata*100,  h , 'replicate');
    h1=figure; imagesc(curvdatafiltered(:,:)); colormap(cmap); colorbar; caxis([-10 10]);
    if illu_status
        hold on; plot([origin-2*radius,origin+worm_length],[j1,j1],'c-');
        hold on; plot([origin-2*radius,origin+worm_length],[j2,j2],'c-');
    end
    origin_y=mcd(istart+j1).IllumRectOrigin(2);
    radius_y=mcd(istart+j1).IllumRectRadius(2);

    if origin_y<=10
        segname=' head';
    elseif origin_y>10 && origin_y<=20
        segname=' neck';
    elseif origin_y>30 && origin_y<=60
        segname=' mid-body';
    elseif origin_y>70
        segname=' tail';
    else
        segname=' various';
    end
    %segname=' head';
    hold on; 
    if illu_status
        if j1
           plot([initiation,initiation,origin_y+radius_y,origin_y+radius_y,initiation],[j1,j2,j2,j1,j1] ,'color',[0.5 0.5 0.5],'linewidth',2);
        end
    end
    title('cuvature diagram');
    set(gca,'XTICK',[1 20 40 60 80 100]);
    set(gca,'XTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);    
    time_auto{cycle,1}=time_auto{cycle,1}-t1;
    %set(gca,'YTICK',1:2*fps:numframes);
    %y_tick=get(gca,'YTICK');
    %set(gca,'YTICKLABEL',time_auto{cycle,1}(y_tick,1));
    set(gca,'YTICK',1:50:numframes);
    set(gca,'YTICKLABEL',istart:50:iend);
    xlabel('fractional distance along the centerline (head=0; tail=1)');
    ylabel('time (s)');

    saveas(h1,strcat(vedio_sequence_s,segname,' cur'),'fig')
    saveas(h1,strcat(vedio_sequence_s,segname,' cur'),'jpg')

    head_curv=zeros(numframes,1);
    medial_curv=zeros(numframes,1);
    tail_curv=zeros(numframes,1);


    %answer = inputdlg({'origin', 'radius'}, '', 1, {num2str(origin), num2str(radius)});
    %origin = str2double(answer{1});
    %radius = str2double(answer{2});

    origin=10;
    radius=8;

    for j=1:numframes
        head_curv(j)=mean(curvdatafiltered(j,origin-radius:origin+radius));
    end
    origin=50;
    radius=10;

    for j=1:numframes
        medial_curv(j)=mean(curvdatafiltered(j,origin-radius:origin+radius));
    end

    origin=80;
    radius=10;

    for j=1:numframes
        tail_curv(j)=mean(curvdatafiltered(j,origin-radius:origin+2*radius));
    end

    segmented_curv=zeros(numframes,5);


    %    head_curv=zeros(numframes,1);
    %	head_curv2=zeros(numframes,1);
    %	head_curv3=zeros(numframes,1); 
    %	head_curv4=zeros(numframes,1);
    %	head_curv5=zeros(numframes,1);
    %	head_curv6=zeros(numframes,1);



    %answer = inputdlg({'origin', 'radius'}, '', 1, {num2str(origin), num2str(radius)});
    %origin = str2double(answer{1});
    %radius = str2double(answer{2});
    %{
    origin=10;
    radius=8;

    segment_length=20;

    start_segment=1:segment_length:100;

    Num_segments=length(start_segment);

    segmented_curv=zeros(numframes,Num_segments);

    for i=1:Num_segments

    for j=1:numframes

    segmented_curv(j,i)=mean(curvdatafiltered(j,start_segment(i):start_segment(i)+segment_length-1));

    end
    end



    %{
    switch origin_y
        case 10
            %segdata_head=cell(20,1);
            segdata_head{vedio_sequence,1}=segmented_curv((j1-100):end,:);
        case 50
            %segdata_medial=cell(20,1);
            segdata_medial{vedio_sequence,1}=segmented_curv((j1-100):end,:);
        case 80
            %segdata_tail=cell(20,1);
            segdata_tail{vedio_sequence,1}=segmented_curv((j1-100):end,:);
    end
    %}

    h1=figure (3);



    %No.1
    onflag=1;
    offflag=1;
     for j=2:numframes
        i=istart+(j-1);
        if ~mcd(i).DLPisOn
            plot([time(j-1),time(j)],[segmented_curv(j-1,1),segmented_curv(j,1)],'k','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'rs');
            hold on;
            if offflag-onflag
            plot([time(j),time(j)],[-10,15],'g-','Linewidth',2);
            offflag=0;
            hold on;
            end

        else
            plot([time(j-1),time(j)],[segmented_curv(j-1,1),segmented_curv(j,1)],'k','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'bs');
           hold on;
           if onflag
           plot([time(j),time(j)],[-10,15],'g-','Linewidth',2);
           onflag=0;
           hold on;
           end
        end
     end

     %No.2

     for j=2:numframes
        i=istart+(j-1);
        if ~mcd(i).DLPisOn
            plot([time(j-1),time(j)],[segmented_curv(j-1,2),segmented_curv(j,2)],'r','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'rs');
            hold on;
        else
            plot([time(j-1),time(j)],[segmented_curv(j-1,2),segmented_curv(j,2)],'r','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'bs');
           hold on;

        end
     end

     %No.3

      for j=2:numframes
        i=istart+(j-1);
        if ~mcd(i).DLPisOn
            plot([time(j-1),time(j)],[segmented_curv(j-1,3),segmented_curv(j,3)],'c','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'rs');
            hold on;
        else
            plot([time(j-1),time(j)],[segmented_curv(j-1,3),segmented_curv(j,3)],'c','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'bs');
           hold on;

        end
     end

     %No.4

     for j=2:numframes
        i=istart+(j-1);
        if ~mcd(i).DLPisOn
            plot([time(j-1),time(j)],[segmented_curv(j-1,4),segmented_curv(j,4)],'m','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'rs');
            hold on;
        else
            plot([time(j-1),time(j)],[segmented_curv(j-1,4),segmented_curv(j,4)],'m','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'bs');
           hold on;

        end
     end

     %No.5

      for j=2:numframes
        i=istart+(j-1);
        if ~mcd(i).DLPisOn
            plot([time(j-1),time(j)],[segmented_curv(j-1,5),segmented_curv(j,5)],'b','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'rs');
            hold on;
        else
            plot([time(j-1),time(j)],[segmented_curv(j-1,5),segmented_curv(j,5)],'b','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'bs');
           hold on;

        end
      end
    title('body curvature diagram')
    xlabel('time(s)')
    ylabel('curvature')
    saveas(h1,strcat(vedio_sequence_s,segname,' seg_cur'),'fig')
    saveas(h1,strcat(vedio_sequence_s,segname,' seg_cur'),'jpg')
    %legend('0-20','20-40','40-60','60-80','80-100');
    %}

      %{
      figure(4)
    for j=2:numframes
        i=istart+(j-1);
        if ~mcd(i).DLPisOn
            plot([time(j-1),time(j)],[head_curv(j-1),head_curv(j-1)],'g','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'rs');
            hold on;
        else
            plot([time(j-1),time(j)],[head_curv(j-1),head_curv(j-1)],'g','Linewidth',2);
            %plot(j,curv_dot/abs(max(curv_dot)),'bs');
           hold on;

        end
     end
      %}
    if j1~=0

        if j2~=0
            frames=length(j1:j2);

            %disp([std(head_curv(max(1,j1-frames):j1)), std(head_curv(j1:j2)),std(head_curv(j2:min(j2+frames,iend-istart+1))),std(head_curv(j2+frames:min(j2+2*frames,iend-istart+1))),std(head_curv(j2:min(j2+2*frames,iend-istart+1)))]);

            disp('head curvature before illu')
            disp(std(head_curv(max(1,j1-3*frames):max(1,j1-2*frames))));

            disp(std(head_curv(max(1,j1-2*frames):max(1,j1-frames))));

            disp(std(head_curv(max(1,j1-frames):j1)));
            disp('medail curvature before illu')
            disp(std(medial_curv(max(1,j1-3*frames):max(1,j1-2*frames))));

            disp(std(medial_curv(max(1,j1-2*frames):max(1,j1-frames))));

            disp(std(medial_curv(max(1,j1-frames):j1)));
            disp('tail curvature before illu')
            disp(std(tail_curv(max(1,j1-3*frames):max(1,j1-2*frames))));

            disp(std(tail_curv(max(1,j1-2*frames):max(1,j1-frames))));

            disp(std(tail_curv(max(1,j1-frames):j1)));
            disp('head curvature during illu');
            disp(std(head_curv(j1:j2)));
            disp('medial curvature during illu');
            disp(std(medial_curv(j1:j2)));
            disp('tail curvature during illu');
            disp(std(head_curv(j1:j2)));
            disp('head curvature after illu');
            disp(std(head_curv(j2:min(j2+frames,iend-istart+1))));
            disp('medial curvature after illu');
            disp(std(medial_curv(j2:min(j2+frames,iend-istart+1))));
            disp('tail curvature after illu');
            disp(std(tail_curv(j2:min(j2+frames,iend-istart+1))));



            disp('start frame and end frame');

            disp([istart iend]);

        else

            N=length(head_curv);
            frames=N-j1+1;
            disp([std(head_curv(max(1,j1-frames):j1)), std(head_curv(j1:end))]);
            disp([istart iend]);
        end

    else 
        disp(std(head_curv));
        disp([istart iend]);
    end


    body_amp1=zeros(1,100);
    body_amp2=zeros(1,100);
    body_amp3=zeros(1,100);
    %% calculate and plot body amplitude
    if illu_status
        for bodyseg=1:100
            if j1==0 && j2==0
                body_amp1(1,bodyseg)=std(curvdatafiltered(:,bodyseg))*2;
            elseif j1~=0 && j2==0
                body_amp1(1,bodyseg)=std(curvdatafiltered(1:j-1,bodyseg))*2;
                body_amp2(1,bodyseg)=std(curvdatafiltered(j1:end,bodyseg))*2;
            elseif j1==0 && j2~=0
                body_amp2(1,bodyseg)=std(curvdatafiltered(1:j2,bodyseg))*2;
                body_amp3(1,bodyseg)=std(curvdatafiltered(j2+1:end,bodyseg))*2;
            elseif j1~=0 && j2~=0
                body_amp1(1,bodyseg)=std(curvdatafiltered(1:j1-1,bodyseg))*2;
                body_amp2(1,bodyseg)=std(curvdatafiltered(j1:j2,bodyseg))*2;
                body_amp3(1,bodyseg)=std(curvdatafiltered(j2+1:end,bodyseg))*2;
            else 
                break
            end
        end
    else
        for bodyseg=1:100
            body_amp1(1,bodyseg)=std(curvdatafiltered(:,bodyseg))*2;
        end
    end
    h1=figure(4);
    plot(1:100,body_amp1,'k',1:100,body_amp2,'g',1:100,body_amp3,'r')
    title('body amplitude diagram')
    xlabel('body coordinate')
    ylabel('amplitude')
    legend('before','during','after');
    saveas(h1,strcat(vedio_sequence_s,segname,' amp'),'fig')
    saveas(h1,strcat(vedio_sequence_s,segname,' amp'),'jpg')
    if illu_status
        body_amp(1,:)=body_amp1;
        body_amp(2,:)=body_amp2;
        body_amp(3,:)=body_amp3;
    else
        body_amp=body_amp1;
    end
    % save different targetting regions amplitude data seperately
    switch segname
        case ' head'
            %segdata_head=cell(20,1);
            amp_head{vedio_sequence,1}=body_amp;
        case ' neck'
            %segdata_medial=cell(20,1);
            amp_neck{vedio_sequence,1}=body_amp;
        case ' mid-body'
            %segdata_tail=cell(20,1);
            amp_midbody{vedio_sequence,1}=body_amp;
        case ' tail'
            amp_tail{vedio_sequence,1}=body_amp;
    end

    %worm_data=struct('Start_frame',istart,'End_frame',iend,'Time',time,'Worm_curvature',curvdatafiltered,'center_of_illu',origin,'radius_of_illu',radius,'mean_curvature',mean_curv,'Video_name',vid); 
    t_w_curvdatafiltered{vedio_sequence,1}=curvdatafiltered;
    t_w_curdata{vedio_sequence,1}=curvdata;
    t_w_amp{vedio_sequence,1}=body_amp;
    t_w_angledata{vedio_sequence,1}=angle_data;
    %{
    t_w_aligning(vedio_sequence,1)=vedio_sequence;
    t_w_aligning(vedio_sequence,2)=istart;
    t_w_aligning(vedio_sequence,3)=iend;
    if exist('DLPon_frame','var')
        t_w_aligning(vedio_sequence,4)=DLPon_frame;
        t_w_aligning(vedio_sequence,6)=initiation;
        t_w_aligning(vedio_sequence,7)=origin_y+radius_y;
    else 
        t_w_aligning(vedio_sequence,4)=0;
    end
    if exist('DLPoff_frame','var')
        t_w_aligning(vedio_sequence,5)=DLPoff_frame;
        t_w_aligning(vedio_sequence,6)=initiation;
        t_w_aligning(vedio_sequence,7)=origin_y+radius_y;
    else
        t_w_aligning(vedio_sequence,5)=0;
    end
    disp('the sequence of the worm is ');
    disp(vedio_sequence);
    %}
    clear j1 j2
end
evalc([wormname,'_amp=t_w_amp']);
evalc([wormname,'_curdata=t_w_curdata']);
evalc([wormname,'_curvdatafiltered=t_w_curvdatafiltered']);
evalc([wormname,'_illu=t_w_illu_time']);
clear ans
disp('+++++++++++++++++++');
disp('       DONE        ');
disp('+++++++++++++++++++');
close all
nisi;






