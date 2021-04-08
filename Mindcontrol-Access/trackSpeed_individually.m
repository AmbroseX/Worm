addpath(pwd);
close all;

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
    
    
    
    %figure (1);
   % plot(centerline(1,:),centerline(2,:),'k-');
   % hold on; plot(Head_position(1),Head_position(2),'ro');
   % hold on; plot(Tail_position(1),Tail_position(2),'bo');
	
    axis off; axis equal; hold on;
    df = diff(centerline,1,2); 
    t = cumsum([0, sqrt([1 1]*(df.*df))]); 
    worm_length=worm_length+t(end);
    cv = csaps(t,centerline,spline_p);
    
   %figure(1);
   % fnplt(cv, '-g'); hold off;   
    
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


worm_length=worm_length/numframes;
%% starts speed tracking





vidfile=strcat(fname(1:end-5),'_HUDS','.avi');
%istart=input(' which frame do you want to start tracking');
%iend=input(' and to which frame to end');
stagepos=trackShmutz(vidfile,istart,iend);


videoDecimationConversion=2;

%Choose a point on the worm to track. Perhapse a point in the neck?
ptAlongWorm=20; %15 percent alont the worm's body length
XX=1;
YY=2;

N=numframes;


com=[mean(Centerline(:,:,XX),2), mean(Centerline(:,:,YY),2)];




%pos=videoDecimationConversion.*[-stagepos(1:N,1), stagepos(1:N,2)]...
%   +[Centerline(1:N,ptAlongWorm,XX), -Centerline(1:N,ptAlongWorm,YY)];

  pos=videoDecimationConversion.*[-stagepos(1:N,1), stagepos(1:N,2)]...
      +[com(1:N,XX), -com(1:N,YY)];

M=30;
intpos=pos;
%intpos=lowpass1d(pos,M)./scaleFactor; %run through a low pass filter
intpos=[intpos(:,1)-intpos(1,1), intpos(:,2)-intpos(1,2)]; %subtract off zero coordinate
diffpos=pos(M+1:end,:)-pos(1:end-M,:);

diffdist=sqrt(sum(diffpos.^2,2));

difftime=time(M+1:end,:)-time(1:end-M,:);


actual_speed=diffdist./difftime/worm_length;




figure;
hold on;
lw=4;
plot(intpos(:,1),intpos(:,2),'bo');
plot(intpos(1,1),intpos(1,2),'ok','LineWidth',lw);

disp('mean actual speed of the worm (L/sec)');
disp(mean(actual_speed));
fprintf('\nspeed tracking period is between %d to %d\n',istart,iend);
%quiver(intpos(1:M./2:end-M,1),intpos(1:M./2:end-M,2),diffpos(1:M./2:end,1),diffpos(1:M/2:end,2),'b');
