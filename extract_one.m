workpath=fullfile('G:','Data','WenLab','Worm_Embed');
addpath(genpath(fullfile(workpath,'libwen')));
filepath='testDLP';
pathname=fullfile(workpath,'rawdata',filepath); %the rawdata's path
yamlfiles = dir(fullfile(pathname,'*.yaml'));


% start_yaml = input('From which file start # of *.yaml file:'); %you can input the star file number
% if isempty(start_yaml)  %if the input is empty then start from the first file
%     start_yaml = 1;
% end
s_yaml=3;

filename = yamlfiles(s_yaml).name;
fname=fullfile(pathname,filename);   % the full pathe of the *.yaml
namepattern = 'w\d*\w*\.yaml';
timepattern = '\d*_\d\d\d\d_';
shortname = regexp(filename,namepattern,'match');
wormname = shortname{1}(1:end-5);



mcd = Mcd_Frame;
mcd = mcd.yaml2matlab(fname);    % a=mcd(x)


frames_afterDLPon =300;
numcurvpts = 100;
proximity = 50;
spline_p = 0.0005;
flip=0;

%delete DLP=0
framnum=length(mcd);
framnum2=0;
for i=1:framnum
    dlp=mcd(i).DLPisOn;
    if dlp==0
        continue;
    elseif dlp==1
        if i<=frames_afterDLPon
            continue;
        elseif i>frames_afterDLPon && i+frames_afterDLPon <framnum
            if mcd(i-frames_afterDLPon).DLPisOn == 1 && mcd(i+frames_afterDLPon).DLPisOn ==1
                framnum2=framnum2+1;
                mcd2(framnum2)=mcd(i);
            end
        else
            continue;
        end
    end  
end

mcd=mcd2;framnum=length(mcd);

clear mcd2

%plot the frame
% for i=1:framnum
%    x(i)=mcd(i).FrameNumber; 
% end
% plot([1:framnum],x)


%cal curve,angle
wormdata.name=yamlfiles(s_yaml).name;
wormdata.wormname=wormname;
wormdata.curv_data=zeros(framnum,numcurvpts);
wormdata.angle_data=zeros(framnum,numcurvpts+1);
wormdata.time=zeros(framnum,1);
wormdata.Centerline=zeros(framnum,100,2);
Head_position=mcd(1).Head;
Tail_position=mcd(1).Tail;
worm_length=0;  %body length in terms of pixels

t1=0;j1=0; j2=0;
for i=1:framnum
    if (norm(mcd(i).Head-Head_position)> norm(mcd(i).Tail-Head_position)) %%head and tail flips
        if norm(mcd(i).Head-Tail_position)<=proximity && norm(mcd(i).Tail-Head_position)<=proximity  %%if the tip points are identified
            flips=~flip;
            Head_position=mcd(i).Tail;
            Tail_position=mcd(i).Head;
        end
    else
        flips = flip;
        Head_position = mcd(i).Head;
        Tail_position = mcd(i).Tail;
    end
    
    if norm(mcd(i).Head-mcd(i).Tail)>proximity
        centerline=reshape(mcd(i).SegmentedCenterline,2,[]);
        if flips
            centerline(1,:)=centerline(1,end:-1:1);
            centerline(2,:)=centerline(2,end:-1:1);
        end
    end
    wormdata.Centerline(i,:,1)=centerline(1,:);
    wormdata.Centerline(i,:,2)=centerline(2,:);
    wormdata.time(i)=mcd(i).TimeElapsed;
    
    df = diff(centerline,1,2);
    t = cumsum([0, sqrt([1 1]*(df.*df))]);%here [0,[1:100]] adds one column by the head, thus the matrix becomes [0:101]
    worm_length=worm_length+t(end);
    cv = csaps(t,centerline,spline_p);
    
    cv2 =  fnval(cv, t)';
    df2 = diff(cv2,1,1); df2p = df2';
    
    splen = cumsum([0, sqrt([1 1]*(df2p.*df2p))]);
    cv2i = interp1(splen+.00001*[0:length(splen)-1],cv2, [0:(splen(end)-1)/(numcurvpts+1):(splen(end)-1)]);
    
    df2 = diff(cv2i,1,1);
    atdf2 =  unwrap(atan2(-df2(:,2), df2(:,1)));
    wormdata.angle_data(i,:) = atdf2';
    
    curv = unwrap(diff(atdf2,1));
    wormdata.curv_data(i,:) = curv';
    
end



clearvars -except wormdata workpath filepath filename yamlfiles
savename=strrep(filename,'.yaml','.mat');
savefolder=fullfile(workpath,'data',filepath);
if exist(savefolder)==0
    mkdir(savefolder);
else
    disp('dir is exist');
end
save(fullfile(savefolder,savename),'wormdata')








