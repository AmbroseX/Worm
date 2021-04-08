addpath(pwd);
close all;
segname='_';

if exist('pathname', 'var')
        try
            if isdir(pathname)
            cd(pathname);
            end
        end
 end
 [filename,pathname]  = uigetfile({'*.yaml'});  
  fname = [pathname filename];
 
 if ~exist('mcd','var')
     mcd=Mcd_Frame;
     mcd=mcd.yaml2matlab(fname);
 end
fprintf('\n%s\n',filename);



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
vedio_part=input('number of vedio part\n','s');

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

origin_y=mcd(istart+j1).IllumRectOrigin(2);
radius_y=mcd(istart+j1).IllumRectRadius(2);

hold on; 
if j2
   plot([origin_y-radius_y,origin_y-radius_y,origin_y+radius_y,origin_y+radius_y,origin_y-radius_y],[j1,j2,j2,j1,j1] ,'color',[0.5 0.5 0.5],'linewidth',2);
end

title('cuvature diagram');


set(gca,'XTICK',[1 20 40 60 80 100]);
set(gca,'XTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);

t3=mcd(istart).TimeElapsed;
time=time-t3;

%set(gca,'YTICK',1:2*fps:numframes);
y_tick=get(gca,'YTICK');
set(gca,'YTICKLABEL',time(y_tick));

xlabel('fractional distance along the centerline (head=0; tail=1)');
ylabel('time (s)');

saveas(h1,strcat(vedio_part,segname,' cur'),'fig')
saveas(h1,strcat(vedio_part,segname,' cur'),'jpg')

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

sequence_1=input('input the sequence of segment data');


 segdata{sequence_1,1}=segmented_curv;



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
saveas(h1,strcat(vedio_part,segname,' seg_cur'),'fig')
saveas(h1,strcat(vedio_part,segname,' seg_cur'),'jpg')
%legend('0-20','20-40','40-60','60-80','80-100');

  
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
        
        disp('head curvature before illumination')
        disp(std(head_curv(max(1,j1-3*frames):max(1,j1-2*frames))));
		
        disp(std(head_curv(max(1,j1-2*frames):max(1,j1-frames))));
        
        disp(std(head_curv(max(1,j1-frames):j1)));
        disp('medail curvature before illumination')
        disp(std(medial_curv(max(1,j1-3*frames):max(1,j1-2*frames))));
		
        disp(std(medial_curv(max(1,j1-2*frames):max(1,j1-frames))));
        
        disp(std(medial_curv(max(1,j1-frames):j1)));
        disp('tail curvature before illumination')
        disp(std(tail_curv(max(1,j1-3*frames):max(1,j1-2*frames))));
		
        disp(std(tail_curv(max(1,j1-2*frames):max(1,j1-frames))));
        
        disp(std(tail_curv(max(1,j1-frames):j1)));
        disp('head curvature during illumination');
        disp(std(head_curv(j1:j2)));
        disp('medial curvature during illumination');
        disp(std(medial_curv(j1:j2)));
        disp('tail curvature during illumination');
        disp(std(head_curv(j1:j2)));
        disp('head curvature after illumination');
        disp(std(head_curv(j2:min(j2+frames,iend-istart+1))));
        disp('medial curvature after illumination');
        disp(std(medial_curv(j2:min(j2+frames,iend-istart+1))));
        disp('tail curvature after illumination');
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


body_amp1=zeros(1,99);
body_cur1=zeros(j1,99);
body_amp2=zeros(1,99);
body_cur2=zeros(j2-j1+1,99);
body_amp3=zeros(1,99);
body_cur3=zeros(numframes-j2+1,99);

%% calculate and plot body amplitude

for bodyseg=2:99
    
    origin=bodyseg;
    radius=1;
    
    for j=1:j
        body_cur1(j,bodyseg)=mean(curvdatafiltered(j,origin-radius:origin+radius));
    end
    
    
    
end

    for bodyseg=1:99
        body_amp1(1,bodyseg)=std(body_cur1(:,bodyseg))*2;
       
        
    end
    

h1=figure(4);


plot(1:99,body_amp1,'k')
title('body amplitude diagram')
xlabel('body coordinate')
ylabel('amplitude')

saveas(h1,strcat(vedio_part,segname,' amp'),'fig')
saveas(h1,strcat(vedio_part,segname,' amp'),'jpg')
body_amp(1,:)=body_amp1;




%worm_data=struct('Start_frame',istart,'End_frame',iend,'Time',time,'Worm_curvature',curvdatafiltered,'center_of_illumination',origin,'radius_of_illumination',radius,'mean_curvature',mean_curv,'Video_name',vid); 

amp{sequence_1,1}=body_amp;







