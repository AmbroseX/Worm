addpath(pwd);

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



if exist('istart', 'var')
        answer = inputdlg({'Start frame', 'End frame'}, 'Cancel to clear previous', 1, ...
            {num2str(istart),num2str(iend)});
    else
        answer = inputdlg({'Start frame', 'End frame'}, '', 1);
end
    
if isempty(answer)
    answer = inputdlg({'Start frame', 'End frame'}, '', 1);
end
    
istart = str2num(answer{1}); %start frame
iend = str2num(answer{2});   %end frame
numframes=iend-istart+1;
numcurvpts=100;
k_buffer=zeros(1,10);  %a dynamic queue for storing the mean curvature of the head
curv=zeros(1,numframes);
curv_dot=zeros(1,numframes);

time=zeros(numframes,1);



for j=1:numframes
    
    i = istart + (j - 1);
    
    centerline=reshape(mcd(i).SegmentedCenterline,2,[]);
    
    time(j)=mcd(i).TimeElapsed;
    
    if mcd(i).DLPisOn && ~mcd(i-1).DLPisOn
        t0=time(j);
    end
    
    
    figure (1);
    plot(centerline(1,:),centerline(2,:),'k-');
    axis equal; hold on;
        
   
    centerline2=centerline;
    
       
    centerline2(2,:)=smooth(centerline(2,:),20);
    centerline2(1,:)=smooth(centerline(1,:),20); %smooth the centerline using the simplist running average algorithm in MATLAB
    
    hold on; plot(centerline2(1,:),centerline2(2,:)); hold off;   
    
    
    [k,k_buffer]=extract_curvature(centerline2(:, 5:30)',k_buffer); %calculate the mean curvature of the head
    k_dot=diff_curv(k_buffer); %calculate the derivative of curvature
    curv(j)=k*numcurvpts;
    curv_dot(j)=k_dot*numcurvpts;
    
end

figure (2);

time=time-t0;

for j=2:numframes
    i=istart+(j-1);
    if ~mcd(i).DLPisOn
        plot([time(j-1),time(j)],[curv(j-1),curv(j)],'k','Linewidth',2);
        %plot(j,curv_dot/abs(max(curv_dot)),'rs');
        hold on;
    else
        plot([time(j-1),time(j)],[curv(j-1),curv(j)],'b','Linewidth',2);
        %plot(j,curv_dot/abs(max(curv_dot)),'bs');
        hold on;
    end
end
        








