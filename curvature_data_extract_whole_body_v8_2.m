tic
%目的：将rawdata中一个文件夹下所有*.yaml文件用DLP的标记，提取出name,头尾位置，angle_data,curv_data
%然后将其导入到data对应文件夹下面
%to get the angle_data and curve_data one by one
%for Rongkang desktop-3070  & Laptap
addpath(genpath(fullfile('G:','Data','WenLab','Worm_Embed','libwen')));
pathname=fullfile('G:','Data','WenLab','Worm_Embed','rawdata','testDLP');
close all;



%for the server-2080Ti
%addpath(genpath(fullfile('/','home','wenlab','xrk','Worm_Embed','libwen')))
%pathname=fullfile('/','home','wenlab','xrk','Worm_Embed','rawdata')
%to test if the pathname existed in the workspace var.
yamlfiles = dir(fullfile(pathname,'*.yaml'));  % get all the *.yaml file's info:name folder
start_yaml = input('From which file start # of *.yaml file:'); %you can input the star file number


if exist('pathname', 'var')
    try
        if isfolder(pathname)
            cd(pathname);
        end
    end
end



if ~exist('mcd','var')
    mcd=Mcd_Frame;
    mcd=mcd.yaml2matlab(fname);
end


%if there is 'istart' in var-space then change it using dialog.
if exist('istart', 'var')
    answer = inputdlg({'Start frame', 'End frame', 'spline fit parameter','flip head and tail?'}, 'Cancel to clear previous', 1, ...
        {num2str(istart),num2str(iend),num2str(spline_p),num2str(flip)});
else
    answer = inputdlg({'Start frame', 'End frame','spline fit parameter','flip head and tail?'}, '', 1,{'0','27000','0.0005','0'});
end

if isempty(answer)
    answer = inputdlg({'Start frame', 'End frame','spline fit parameter','flip head and tail?'}, '', 1,{'0','27000','0.0005','0'});
end

istart = str2num(answer{1});
iend = str2num(answer{2});
spline_p = str2num(answer{3});
flip = str2num(answer{4});  % if the head and tail do not recognize the opposite,flip=0,otherwise flip=1

numframes=iend-istart+1; %the totle number will analysis

numcurvpts=100;

proximity = 50;

curvdata=zeros(numframes,numcurvpts);
angle_data = zeros(numframes,numcurvpts+1);
time=zeros(numframes,1);

Head_position=mcd(istart).Head;
Tail_position=mcd(istart).Tail;

worm_length=0;  %body length in terms of pixels

mcd2=mcd;

t1=0;

j1=0; j2=0;

Centerline=zeros(numframes,100,2);
clear info;
for j=1:numframes
    
    i = istart + (j - 1);
    
    if (norm(mcd(i).Head-Head_position)> norm(mcd(i).Tail-Head_position)) %%head and tail flips
        if norm(mcd(i).Head-Tail_position)<=proximity && norm(mcd(i).Tail-Head_position)<=proximity  %%if the tip points are identified
            flip=~str2num(answer{4});
            Head_position=mcd(i).Tail;
            Tail_position=mcd(i).Head;
            %mcd2(i).Head=Head_position;
            %mcd2(i).Tail=Tail_position;
        end
    else
        flip = str2num(answer{4});
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
    
    
    
    
    time(j)=mcd(i).TimeElapsed;
    %{
    if mcd(i).DLPisOn && ~mcd(i-1).DLPisOn
        t1=time(j);
		w1=t1;
        j1=j;
        %origin=100-mcd(i).IllumRectOrigin(2);
        %radius=mcd(i).IllumRectRadius(2);
    end
    if ~mcd(i).DLPisOn && mcd(i-1).DLPisOn
        t2=time(j);
        w2=t2;
		j2=j;
    end
    %}
    
    
    %figure (1);
    %plot(centerline(1,:),centerline(2,:),'k-');
    %hold on; plot(Head_position(1),Head_position(2),'ro');
    %hold on; plot(Tail_position(1),Tail_position(2),'bo');
    
    %axis off; axis equal; hold on;
    %calculate the length of worm using the vector coordinates of each point of the worm
    df = diff(centerline,1,2);
    t = cumsum([0, sqrt([1 1]*(df.*df))]);%here [0,[1:100]] adds one column by the head, thus the matrix becomes [0:101]
    worm_length=worm_length+t(end);
    cv = csaps(t,centerline,spline_p);
    
    
    
    
    %figure(1);
    %fnplt(cv, '-g'); hold off;
    
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

%use the wormlength of the first frame as a constent, then when
%head-tail distence is lower than a special vlaue, we then consider it
%an omega turn



for m=1:numframes
    n=istart+m-1;
    centerline_new=reshape(mcd(n).SegmentedCenterline,2,[]);
    if flip
        centerline(1,:)=centerline(1,end:-1:1);
        centerline(2,:)=centerline(2,end:-1:1);
    end
    dif=diff(centerline_new,1,2);
    
    sumdif=cumsum(sqrt([1 1]*(dif.*dif)));
    wlength=sumdif(end);
    
    distance_h2t=norm(centerline_new(:,2)-centerline_new(:,end));
    if wlength<distance_h2t
        bool=0;
    else
        bool=1;
    end
    info(m,1)=n;
    info(m,2)=distance_h2t;
    info(m,3)=wlength;
    info(m,4)=bool;
    
    
end

%{
cmap=redgreencmap;
cmap(:,3)=cmap(:,2);
cmap(:,2)=0;
origin=10;
radius=8;
%}
worm_length=worm_length/numframes;

answer = inputdlg({'time filter', 'body coord filter', 'mean=0, median=1'}, '', 1, {num2str(5), num2str(10), '0'});
timefilter = str2double(answer{1});
bodyfilter = str2double(answer{2});

h = fspecial('average', [timefilter bodyfilter]);
curvdatafiltered = imfilter(curvdata*100,  h , 'replicate');


figure;
imagesc(curvdatafiltered(:,:)); colormap(map); colorbar; caxis([-10 10]);
hold on; plot([origin-2*radius,origin+worm_length],[j1,j1],'c-');
hold off
figure
hold on; plot([origin-2*radius,origin+worm_length],[j2,j2],'c-');
hold off
%hold on; plot([origin-radius,origin-radius,origin+radius,origin+radius,origin-radius],[j1,j2,j2,j1,j1] ,'color',[0.5 0.5 0.5],'linewidth',2);
title('cuvature diagram');
set(gca,'XTICK',[1 20 40 60 80 100]);
set(gca,'XTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);
time=time-t1;
%set(gca,'YTICK',1:2*fps:numframes);
y_tick=get(gca,'YTICK');
set(gca,'YTICKLABEL',time(y_tick));
xlabel('fractional distance along the centerline (head=0; tail=1)');
ylabel('time (s)');
head_curv=zeros(numframes,1);
%answer = inputdlg({'origin', 'radius'}, '', 1, {num2str(origin), num2str(radius)});
%origin = str2double(answer{1});
%radius = str2double(answer{2});
origin=10;
radius=8;

for j=1:numframes
    head_curv(j)=mean(curvdatafiltered(j,origin-radius:origin+radius));
end

for j=1:numframes
    head_curv(j)=mean(curvdatafiltered(j,origin+radius:origin+2*radius));
end
figure (3);
for j=2:numframes
    i=istart+(j-1);
    if ~mcd(i).DLPisOn
        plot([time(j-1),time(j)],[head_curv(j-1),head_curv(j)],'k','Linewidth',2);
        %plot(j,curv_dot/abs(max(curv_dot)),'rs');
        hold on;
    else
        plot([time(j-1),time(j)],[head_curv(j-1),head_curv(j)],'g','Linewidth',2);
        %plot(j,curv_dot/abs(max(curv_dot)),'bs');
        hold on;
    end
end
if j1~=0
    if j2~=0
        frames=length(j1:j2);
        %disp([std(head_curv(max(1,j1-frames):j1)), std(head_curv(j1:j2)),std(head_curv(j2:min(j2+frames,iend-istart+1))),std(head_curv(j2+frames:min(j2+2*frames,iend-istart+1))),std(head_curv(j2:min(j2+2*frames,iend-istart+1)))]);
        disp('head curvature before illumination')
        disp(std(head_curv(max(1,j1-3*frames):max(1,j1-2*frames))));
        disp(std(head_curv(max(1,j1-2*frames):max(1,j1-frames))));
        disp(std(head_curv(max(1,j1-frames):j1)));
        disp('head curvature during illumination');
        disp(std(head_curv(j1:j2)));
        disp('head curvature after illumination');
        disp(std(head_curv(j2:min(j2+frames,iend-istart+1))));
        disp(std(head_curv(j2+frames:min(j2+2*frames,iend-istart+1))));
        disp(std(head_curv(j2+2*frames:min(j2+3*frames,iend-istart+1))));
        disp(std(head_curv(j2+3*frames:min(j2+4*frames,iend-istart+1))));
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
worm_data=struct('Start_frame',istart,'End_frame',iend,'Time',time,'Worm_curvature',curvdatafiltered,'center_of_illumination',origin,'radius_of_illumination',radius,'mean_curvature',mean_curv,'Video_name',vid);
