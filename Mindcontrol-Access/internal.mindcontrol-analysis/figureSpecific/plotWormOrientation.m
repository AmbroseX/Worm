function plotWormOrientation(file,m)
%plotWormOrientation(file,m)
%
%Load a file exported from the YAML previewer software and display the
% worm's average orientation as a vector over time.
%
% This takes the average of the orientation of all of the worm's body
% segments. So if exactly half of the worm's body segments are pointed
% North and half are pointed south, the resulting vector would be 0. 
%
% The magnitude says something about the portion of the worm oriented in
% that direction.
%
% The vector is normalized by the worm's instantaneous arclength
%
% m specifices only display every mth quiver. m=10 by default
%
% Note this requries the filledRibbon function that Marc Gershow wrote:
% http://stackoverflow.com/questions/3979172/plotting-evolution-of-2d-vecto
% r-in-3d-as-a-ribbon-in-matlab/3982099#3982099
%
% Andrew Leifer
% leifer@fas.harvard.edu
% 18 October 2010


%Load a file that was exported from the previewer software
load(file);

if ~exist('filledRibbon','file')
   msgbox('filledRibbon could not be found. Please add it to the path or download it from http://stackoverflow.com/questions/3979172/plotting-evolution-of-2d-vector-in-3d-as-a-ribbon-in-matlab/3982099#3982099','Error','error');
   error('The function filledRibbon() could not be found. Please add it to the path or download it from http://stackoverflow.com/questions/3979172/plotting-evolution-of-2d-vector-in-3d-as-a-ribbon-in-matlab/3982099#3982099');
end

% copy the x and y coordinates of the worm's centerline (recall 0,0 is
% head)
xpos=handles.SegmentedCenterlinex_data;
ypos=handles.SegmentedCenterliney_data;

N=size(xpos,1); %number of frames


figure;
hold on;
for k=1:N;
    %take the difference of neighboring points to get tangent vectors
    xdiff=diff(xpos(k,:));
    ydiff=diff(ypos(k,:));
    
    %For each tangent vector sum the individual vector lengths
    instWormlength= sum((xdiff.^2+ydiff.^2).^0.5);
    
    %sum the individual components to get the overall heading of the worm
    a(1,k)=sum(xdiff)/instWormlength;
    a(2,k)=sum(ydiff)/instWormlength;
    quiver(a(1,k),a(2,k));
    
end
x=0:N-1;
y=zeros([1,N]);
z=zeros([1,N]);
u=zeros([1,N]);
v=a(1,:);
w=a(2,:);
c='r';
figure;
lim=max(max(a));
ylim([-lim,lim]);
zlim([-lim,lim]);

% Marc Gershow's ribbon plotting function
h=filledRibbon (x,y,z,u,v,w,c, 'FaceAlpha',0.2);
gcf;
hold on;

if ~exist('m','var')
    m=10; %only display every mth quiver
end

if isempty(m)
    m=10;
end

%add the quiver (but don't plot allof the quivers because that is just too
%busy
quiver3(x(1:m:end),y(1:m:end),z(1:m:end),u(1:m:end),v(1:m:end),w(1:m:end),0);
axis vis3d;
    