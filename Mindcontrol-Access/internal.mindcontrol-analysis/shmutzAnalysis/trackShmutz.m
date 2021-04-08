function shmutzpos=trackShmutz(vidfile,istart,iend)
% shmutzpos = trackShmutz(vidfile)
%
% Given a video file, this function will track a bright white piece of
% shmutz throughout the vidoe and will return the position of the piece of
% shmutz at each frame.
%
% The algorithm is extermely simple. The user is asked to draw a circular
% region. In each frame the script returns the location of the brigthest
% pixel within the boundary of the circle. At each frame the circle is
% recentered on this point. If more than one pixel has the same brightest
% value, one is chosen arbitrarily.
%
% This seems to work perfectly for to accurately the stage motion in
% videos taken on the CoLBeRT setup with MindControl software.
%
% Wrilltten by Andrew Leifer
% leifer@fas.harvard.edu
% Released under the GNU General Public License

%Read in Video
playstep = 20; % update the ploted figure each playstep interval.
if ~exist(vidfile,'file')
    error('Could not find specified video file!');
end

vid=VideoReader(vidfile);
I=read(vid,istart);
%Get number of frames
N=get(vid,'NumberOfFrames');

if (N<2)
    error('There does not seem to be enough frames in this video to do anythign useful.');
end

mainfig=figure;
h_I=imshow(I);
a_I=get(h_I,'Parent');
roi=imellipse;
roi_pos=getPosition(roi);


N=iend-istart+1;
brightest=zeros(N,2);

for k=1:N
    
%Create a Mask
I=read(vid,k+istart-1);
BW = createMask(roi,h_I);

%Multiply the binary mask with our image
%and find the brightest point
masked=rgb2gray(I).*uint8(BW);
value=max(max(masked));
[r,c]=find(masked==value);


brightest(k,:)=[r(1) c(1)]; %save the value

if mod(k,playstep) == 1
    figure(mainfig)
    %Update the figure with the newest frame
    h_I=imshow(I,'Parent',a_I);

    %Plot this point and all subsequent points
    hold on;
    plot(brightest(1:end,2),brightest(1:end,1),'ro');
    hold off;
end
%Reinstate the ellipse
roi=imellipse(gca,roi_pos);


%get current region
roi_pos=getPosition(roi);

%Set the position of the future ellipse to be the same as the tracked point
setPosition(roi,[brightest(k,2)-roi_pos(3)/2, brightest(k,1)-roi_pos(4)/2, roi_pos(3), roi_pos(4) ]);
end
close;

shmutzpos=[brightest(:,2),brightest(:,1)] ;