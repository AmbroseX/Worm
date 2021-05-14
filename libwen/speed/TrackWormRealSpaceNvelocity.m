function speed = TrackWormRealSpaceNvelocity(segments,start_counting_segment,time_speed,worm_length,com,stagepos)
N = size(stagepos,1);
pos = zeros(N,2,segments-start_counting_segment+1);
for segment = 1:segments-start_counting_segment+1
    pos(:,:,segment) = [-stagepos(1:N,1), stagepos(1:N,2)]...
      +[com(1:N,1,segment), -com(1:N,2,segment)];
end

%% calculate worm velocity
M = 1; % points to be overlapped and smoothed
diffpos = pos(M+1:end,:,:)-pos(1:end-M,:,:);
diffdist=mean(sqrt(sum(diffpos.^2,2)),3);
difftime=time_speed(M+1:end,:)-time_speed(1:end-M,:);

actual_speed=diffdist./difftime/worm_length;

% actual_speed=diffdist./difftime/1000; % change unit of measurement from micronpersecond to mm/s;
speed_zero = zeros(M,1)*NaN;
speed = mean(actual_speed,2);
speed = [speed_zero;speed];

%% trajectory of body point(s) 
%{
intpos=mean(pos,3);
intpos=[intpos(:,1)-intpos(1,1), intpos(:,2)-intpos(1,2)];
figure;
hold on;
lw=4;
plot(intpos(:,1),intpos(:,2),'bo');
plot(intpos(1,1),intpos(1,2),'ok','LineWidth',lw);
hold off
quiver(intpos(1:M/2:end-M,1),intpos(1:M/2:end-M,2),diffpos(1:M/2:end,1),diffpos(1:M/2:end,2),'b');
%}
end