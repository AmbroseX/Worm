clear flip
close all
left_bound = 5;
tremor_thshd = 1;
answers = inputdlg({'trail number', 'region (head = 1, mid = 2, tail = 3','flip number'}, '',1);
trail = str2num(answers{1});
region = str2num(answers{2});
flip_se = str2num(answers{3});
failed = 0; %no whole body effect
successed = 0; %whole body effect
head_stopage = cell(4,length(trail));
%if exist('frames_beforeDLPon','var')
%   DLPonStartPoint = frames_beforeDLPon;
%else
%    DlPon_point = 350;
%end
if ~exist('t_w_illu_time','var')
    t_w_illu_time = this_worm_illumination_time;
end
if ~exist('t_w_curvdatafiltered','var') && exist('this_worm_curvdatafiltered','var')
    t_w_curvdatafiltered = this_worm_curvdatafiltered;
end  

for i = 1:length(trail) % i is the sequence of trail 
    t_w_curvdatafiltered{trail(i),1}(find((t_w_curvdatafiltered{trail(i),1} > -tremor_thshd).*( t_w_curvdatafiltered{trail(i),1} < tremor_thshd))) = 0;
    DLPonStartPoint = find(time_auto{trail(i),1} == 0);
    DLPonEndPoint = DLPonStartPoint+t_w_illu_time(trail(i),2)-t_w_illu_time(trail(i),1);
    if find(flip_se == trail(i))
        c2 = t_w_curvdatafiltered{trail(i),1}(DLPonStartPoint:DLPonEndPoint,:)>0;
        c2 = flip(c2,2);
    else
        c2 = t_w_curvdatafiltered{trail(i),1}(DLPonStartPoint:DLPonEndPoint,:)>0;
    end
    c3 = edge(single(c2),'sobel',0);
    [c4,numlab] = bwlabel(c3);
    figure;imagesc(c4)
    pause(1);
    for j = 1:numlab % j is the number of labeled line
        %% a curve that starts and ends from left_bound smaller than 6 and last longer than 5
        c5 = c4 == j;
        c5(:,left_bound:end) = 0;
        c5 = edge(c5,'sobel',0);
        [c5,overlap_numlab] = bwlabel(c5);
        if overlap_numlab >= 2
            ymean = zeros(overlap_numlab,1);
            for overlap_se = 1:overlap_numlab
                [y,x] = find(c5 == overlap_se);
                ymean(overlap_se) = mean(y);
            end
            y_distance = zeros(overlap_numlab-1,1);
            for overlap_se = 1:overlap_numlab-1
                y_distance(overlap_se) = ymean(overlap_se+1)-ymean(overlap_se);
            end
            y_distance = mean(y_distance);
        else
            y_distance = 0;
        end
        %% begin the main calculation
        if y_distance >= 5 && y_distance <= 100
            [y,x] = find(c4 == j);
            if region == 2
                if max(x) >= 40
                    %figure;imagesc(c2);
                    %if max(y) >=150
                    %('mid-targeting did not function','!! Warning !!')
                    x_max(j) = max(x);
                    failed = failed+1;
                    %end
                else
                    x_max(j) = max(x);
                    successed = successed+1;
                end
            end
       % elseif max(x) <= 40 && (20 <= max(y)-min(y) && max(y)-min(y) <= 100)
        %    x_max(j) = max(x);
            
        end
        
    end
    if exist('x_max','var')
    head_stopage{1,i} = x_max(x_max ~= 0);
    head_stopage{2,i} = mean( head_stopage{1,i});
    head_stopage{3,i} = failed;
    head_stopage{4,i} = successed;
    end
    clear x_max
end
save([wormname,'headstopage'],'head_stopage');
