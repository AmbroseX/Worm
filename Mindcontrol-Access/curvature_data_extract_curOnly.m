addpath(pwd);
close all;

if exist('pathname', 'var')
        try
            if isdir(pathname)
            cd(pathname);
            end
        end
end

%only get ymal file once in one worm's analysis
 if ~exist('fname')
     [filename,pathname]  = uigetfile({'*.yaml'});
 end
 
 fname = [pathname filename];
 
 if ~exist('mcd','var')
     mcd=Mcd_Frame;
     mcd=mcd.yaml2matlab(fname);
 end
fprintf('\n%s\n',filename);

if ~exist('wormname')
    wormname=input('Input the worm name \n','s');
end
% Extrect the DLPon frames of this worm's experiment.
if ~exist('this_worm_illumination_time')
    numframes_total=size(mcd);
    j=1; 
    numframes_total=numframes_total(2);
    this_worm_illumination_time=zeros(30,6);
    numframes_total=size(mcd);  
    numframes_total=numframes_total(2);
    j=1; 
    this_worm_illumination_time=zeros(30,6);
    for i=1:numframes_total-1
         if ~mcd(i).DLPisOn && mcd(i+1).DLPisOn
            origin_y=mcd(i).IllumRectOrigin(2);
            radius_y=mcd(i).IllumRectRadius(2);
            if origin_y-radius_y<0
                 initiation=0;
            else
                initiation=origin_y-radius_y;
            end
            this_worm_illumination_time(j,1)=i;
            this_worm_illumination_time(j,3)=initiation;
            this_worm_illumination_time(j,4)=origin_y+radius_y;
            if i>300
                this_worm_illumination_time(j,5)=i-mod(i,100)-3*100;
            else 
                this_worm_illumination_time(j,5)=i;
            end
         end
        if mcd(i).DLPisOn && ~mcd(i+1).DLPisOn
            this_worm_illumination_time(j,2)=i;
            if i<numframes_total-300
                this_worm_illumination_time(j,6)=i-mod(i,100)+3*100;
            else
                this_worm_illumination_time(j,6)=i;
            end
            j=j+1;
        end  
    end
end

if exist('istart', 'var')
        answer = inputdlg({'Start frame', 'End frame', 'spline fit parameter','flip head and tail?'}, 'Cancel to clear previous', 1, ...
            {num2str(istart),num2str(iend),num2str(spline_p),num2str(flip)});
    else
        answer = inputdlg({'Start frame', 'End frame','spline fit parameter','flip head and tail?'}, '', 1);
end
    
if isempty(answer)
    answer = inputdlg({'Start frame', 'End frame','spline fit parameter','flip head and tail?'}, '', 1);
end
    
istart = str2num(answer{1});
iend = str2num(answer{2});
spline_p = str2num(answer{3});
flip = str2num(answer{4});

numframes=iend-istart+1;

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

%display the target regoin in the screen
origin_y=mcd(istart+j1).IllumRectOrigin(2);
radius_y=mcd(istart+j1).IllumRectRadius(2);
if origin_y-radius_y<0
   initiation=0;
else
    initiation=origin_y-radius_y;
end
for j=1:numframes
    i = istart + (j - 1);
    if ~mcd(i-1).DLPisOn && mcd(i).DLPisOn
        if origin_y-radius_y<0
            initiation=0;
        else
            initiation=origin_y-radius_y;
        end
        fprintf('\n ************************** \n the target region is [%d,%d].\n **************************',initiation,origin_y+radius_y);
        DLPon_frame=i;
    end
    if mcd(i-1).DLPisOn && ~mcd(i).DLPisOn
        if origin_y-radius_y<0
            initiation=0;
        else
            initiation=origin_y-radius_y;
        end
        fprintf('************************** \n the target region is [%d,%d].\n **************************\n',initiation,origin_y+radius_y);
        DLPoff_frame=i;
        break;
    end
end
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
    
    
    
    figure (1);
    plot(centerline(1,:),centerline(2,:),'k-');
    hold on; plot(Head_position(1),Head_position(2),'ro');
    hold on; plot(Tail_position(1),Tail_position(2),'bo');
	
    axis off; axis equal; hold on;
    df = diff(centerline,1,2); 
    t = cumsum([0, sqrt([1 1]*(df.*df))]); 
    worm_length=worm_length+t(end);
    cv = csaps(t,centerline,spline_p);
    
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
vedio_sequence=input('input the sequence of this video part\n');
vedio_sequence_s=num2str(vedio_sequence);
vedio_sequence_s=strcat(wormname,'_',vedio_sequence_s);

cmap=redgreencmap;
cmap(:,3)=cmap(:,2);
cmap(:,2)=0;
origin=10;
radius=8;

worm_length=worm_length/numframes;

answer = inputdlg({'time filter', 'body coord filter', 'mean=0, median=1'}, '', 1, {num2str(5), num2str(10), '0'});
timefilter = str2double(answer{1});
bodyfilter = str2double(answer{2});



h = fspecial('average', [timefilter bodyfilter]);
curvdatafiltered = imfilter(curvdata*100,  h , 'replicate');
h1=figure; imagesc(curvdatafiltered(:,:)); colormap(cmap); colorbar; caxis([-10 10]);

 
hold on; plot([origin-2*radius,origin+worm_length],[j1,j1],'c-');
hold on; plot([origin-2*radius,origin+worm_length],[j2,j2],'c-');


switch origin_y
    case 20
        segname=' neck';
    case 10
        segname=' head';
    case 50
        segname=' medial';
    case 80
        segname=' tail';
end
%segname=' head';
hold on; 
if j2
   plot([origin_y-radius_y,origin_y-radius_y,origin_y+radius_y,origin_y+radius_y,origin_y-radius_y],[j1,j2,j2,j1,j1] ,'color',[0.5 0.5 0.5],'linewidth',2);
end

title('cuvature diagram');


set(gca,'XTICK',[1 20 40 60 80 100]);
set(gca,'XTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);


time=time-t1;

%set(gca,'YTICK',1:2*fps:numframes);
y_tick=get(gca,'YTICK');
set(gca,'YTICKLABEL',time(y_tick));

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
disp([istart iend]);
this_worm_curdata{vedio_sequence,1}=curvdata;
%this_worm_amp{vedio_sequence,1}=body_amp;
this_worm_aligning(vedio_sequence,1)=vedio_sequence;
this_worm_aligning(vedio_sequence,2)=istart;
this_worm_aligning(vedio_sequence,3)=iend;
if exist('DLPon_frame')
    this_worm_aligning(vedio_sequence,4)=DLPon_frame;
    this_worm_aligning(vedio_sequence,6)=initiation;
    this_worm_aligning(vedio_sequence,7)=origin_y+radius_y;
else 
    this_worm_aligning(vedio_sequence,4)=0;
end
if exist('DLPoff_frame')
    this_worm_aligning(vedio_sequence,5)=DLPon_frame;
    this_worm_aligning(vedio_sequence,6)=initiation;
    this_worm_aligning(vedio_sequence,7)=origin_y+radius_y;
else
    this_worm_aligning(vedio_sequence,5)=0;
end
clear DLPon_frame;
clear DLPon_frame;

