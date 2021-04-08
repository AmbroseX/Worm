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
        answer = inputdlg({'Start frame', 'End frame', 'flip head and tail?'}, 'Cancel to clear previous', 1, ...
            {num2str(istart),num2str(iend),num2str(flip)});
    else
        answer = inputdlg({'Start frame', 'End frame','flip head and tail?'}, '', 1);
end
    
if isempty(answer)
    answer = inputdlg({'Start frame', 'End frame','flip head and tail?'}, '', 1);
end
    
istart = str2num(answer{1});
iend = str2num(answer{2});
flip = str2num(answer{3});

numframes=iend-istart+1;

head_tail_distance=zeros(numframes,1);

time=zeros(numframes,1);

k=1;

clear t1;
clear t2;
clear loc;

for j=1:numframes
    
    i = istart + (j - 1);
    
    head_tail_distance(j)=norm(mcd(i).Head-mcd(i).Tail);
    
    time(j)=mcd(i).TimeElapsed;
    
    if mcd(i).DLPisOn && ~mcd(i-1).DLPisOn
        t1(k)=time(j);
        loc(k)=j;
    end
    
    if ~mcd(i).DLPisOn && mcd(i-1).DLPisOn
        t2(k)=time(j);
        j2=j;
        k=k+1;
    end
    
    
end

t=time-time(1);

t1=t1-time(1);
t2=t2-time(1);

r=estimate_turning_rate(head_tail_distance,t);
[counts,binned_t]=count_turns(head_tail_distance,t);

figure; plot(t,r,'k-');
ylim=get(gca,'Ylim');

for k=1:length(t1)
    hold on; plot([t1(k) t1(k)], [ylim(1) ylim(2)],'color','b','linewidth',2);
    hold on; plot([t2(k) t2(k)], [ylim(1) ylim(2)],'color','b','linewidth',2);
end

xlabel('time (s)');
ylabel('turning rate (1/s)');

R=zeros(length(t1),2001);

T=zeros(length(t1),2001);

COUNTS=zeros(length(t1),51);

T_BIN=zeros(length(t1),51);

for k=1:length(t1)
    
    jj=loc(k);
    R(k,:)=r(jj-1000:jj+1000);
    T(k,:)=t(jj-1000:jj+1000)-t1(k);
    [~,loc_bin]=histc(t1(k),binned_t);
    COUNTS(k,:)=counts(loc_bin-25:loc_bin+25);
    T_BIN(k,:)=round((binned_t(loc_bin-25:loc_bin+25)-t1(k)));
    
    
end

figure; plot(mean(T,1),mean(R,1),'k-');

figure; plot(mean(T_BIN,1),mean(COUNTS,1),'k-');

if length(questdlg('Save this data? '))==3
    [fn, savepathname]= uiputfile('*.mat', 'choose file to save', strcat(fname(1:end-5),'.mat'));
    if length(fn) > 1
        fnamemat = strcat(savepathname,fn);
        save(fnamemat);
    end
    
end

    
    
    
    

%worm_data=struct('Start_frame',istart,'End_frame',iend,'Time',time,'Worm_curvature',curvdatafiltered,'center_of_illumination',origin,'radius_of_illumination',radius,'mean_curvature',mean_curv,'Video_name',vid); 








