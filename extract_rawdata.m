tic
%目的：将rawdata中一个文件夹下所有*.yaml文件用DLP的标记，提取出name,头尾位置，angle_data,curv_data
%然后将其导入到data对应文件夹下面
filepathname={'silence','testDLP'};

%for Rongkang desktop-3070  & Laptap
workpath=fullfile('G:','Data','WenLab','Worm_Embed');
%For the 2080Ti
% workpath=fullfile('/','home','wenLab','xrk','Worm_Embed');

addpath(genpath(fullfile(workpath,'libwen')));

for filepath_num=1:length(filepathname)
    filepath=filepathname{filepath_num};
    pathname=fullfile(workpath,'rawdata',filepath); %the rawdata's path
    yamlfiles = dir(fullfile(pathname,'*.yaml'));
    if length(yamlfiles)==0
        disp(strcat(filepath,' folder has no *.yaml'))
        continue;
    end
    disp(strcat('Begain to process:',filepath));
    for s_yaml=1:length(yamlfiles)
        
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
        
        mcd=mcd2;
        framnum=length(mcd);
        
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
        wormdata.TimeElapsed=zeros(framnum,1);
        wormdata.Centerline=zeros(framnum,100,2);
        wormdata.StagePosition=zeros(framnum,2);
        wormdata.StageFeedbackTarget=zeros(framnum,2);
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
            wormdata.TimeElapsed(i)=mcd(i).TimeElapsed;
            wormdata.StagePosition(i,:)=mcd(i).StagePosition;
            wormdata.StageFeedbackTarget(i,:)=mcd(i).StageFeedbackTarget;
            
            df = diff(centerline,1,2); %列差分计算，相邻点做差分
            t = cumsum([0, sqrt([1 1]*(df.*df))]);%求矩阵或向量的累积和，here [0,[1:100]] adds one column by the head, thus the matrix becomes [0:101]
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
        
        clearvars -except wormdata workpath pathname filepathname filepath filename yamlfiles
        savename=strrep(filename,'.yaml','.mat');
        savefolder=fullfile(workpath,'data',filepath);
        
        if exist(savefolder)==0
            disp('dir is not exist');
            mkdir(savefolder);
            disp('make dir success');
        else
            disp('dir is exist');
        end
        
        %calculate the relative speed
        position=wormrelativePosition(wormdata);
        wormdata.speed=wormSpeed(position); % 三列time,Vx,Vy  mm/s
        
        save(fullfile(savefolder,savename),'wormdata')
        disp(['Save file',savename, 'success'])
        
    end
end
toc





