function [speed, stagepos] = queryWormspeed(mcd,worm_length,Stage_Pos,reference_point_mtd,CalculationMethodSelection,varargin)
    
  
    %% if this code is running as a script
    %{
    LoadWormdata;
    if ~exist('running_count','var')
        running_count = 1;
        wormspeed = cell(1);
        worm_stagepos = cell(1);
    end

    if ~exist('fname','var')
         [filename,pathname]  = uigetfile({'*.avi'});
         fname = [pathname filename];
    end
    %}
    %% if this code is running as a function
    if length(varargin) == 2
        istart = varargin{1};
        iend = varargin{2};
        vedioparts = 1;
    elseif length(varargin) == 4
        istart = varargin{1};
        ibreak = varargin{2};
        iend = varargin{3};
        fname = varargin{4};
        vedioparts = 2;
    end
    %% if this code is running as a script
    %{
    checkpoints = input('istart iend or a break point between them\n');
    if length(checkpoints) == 2
        
        iend = checkpoints(2);
       
    elseif length(checkpoints) == 3
        istart = checkpoints(1);
        ibreak = checkpoints(2);
        iend = checkpoints(3);
        
    end
    %}
    micronperunit = 0.05;
    micronperpxl = 2.43; 
    % micronperpxl = 3.16;
    videoDecimationConversion = 2;
    TrackWorm;
    if exist('fname','var')
     if fname(end-3:end) ~= '.'
        fname(end-4:end) = '_HUDS';
        fname = [fname '.avi'];
     end
    end
   
    %reference_point_mtd = 1; % 1 means using stageposition in yaml, 2 means using Shmutz
    switch reference_point_mtd
        case 1
            stagepos = zeros(N,2);
            for i = 1:N
               % stagepos(i,:) = videoDecimationConversion.*(mcd(i+istart-1).StagePosition.*micronperunit/micronperpxl-mcd(i+istart-1).StageFeedbackTarget);
               stagepos(i,:) = videoDecimationConversion.*(Stage_Pos(i+istart-1).*micronperunit/micronperpxl-mcd(i+istart-1).StageFeedbackTarget);
               % stagepos(i,:) = Stage_Pos(i+istart-1).*micronperunit;
            end
            speed = TrackWormRealSpaceNvelocity(segments,start_counting_segment,time_speed,worm_length,com,stagepos);
        case 2
            if vedioparts == 1
                stagepos = videoDecimationConversion.*trackShmutz(fname,istart,iend);
                speed = TrackWormRealSpaceNvelocity(segments,start_counting_segment,time_speed,worm_length,com,stagepos);
            elseif vedioparts == 2
                stagepos1 = videoDecimationConversion.*trackShmutz(fname,istart,ibreak);
                stagepos2 = videoDecimationConversion.*trackShmutz(fname,ibreak-20,iend);
                speed1 = TrackWormRealSpaceNvelocity(segments,start_counting_segment,time_speed(1:ibreak-istart+1),worm_length,com(1:ibreak-istart+1,:,:),stagepos1);
                speed2 = TrackWormRealSpaceNvelocity(segments,start_counting_segment,time_speed(ibreak-istart-19:end),worm_length,com(ibreak-istart-19:end,:,:),stagepos2);
                speed = [speed1;speed2(2:end)];
                stagepos = [stagepos1;stagepos2];
            end
    end
    %% if this code is running as a script
    %{
    wormspeed{running_count} = speed;
    worm_stagepos{running_count} =stagepos;
    running_count = running_count + 1;
    fprintf('%dth analysed\n',running_count-1);
    %}
    
    figure;plot(speed);
end
