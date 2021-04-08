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


%cal curve,angle
wormdata.name=yamlfiles(s_yaml).name;
wormdata.wormname=wormname;
wormdata.curve_data=[];
wormdata.angle_data=[];

%plot the frame
% for i=1:framnum
%    x(i)=mcd(i).FrameNumber; 
% end
% plot([1:framnum],x)










