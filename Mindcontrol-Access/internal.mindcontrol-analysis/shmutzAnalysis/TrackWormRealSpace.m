

videoDecimationConversion=2;

%Choose a point on the worm to track. Perhapse a point in the neck?
ptAlongWorm=20; %15 percent along the worm's body length
XX=1;
YY=2;

N=iend-istart+1;


com=[mean(Centerline(:,:,XX),2), mean(Centerline(:,:,YY),2)];




%pos=videoDecimationConversion.*[-stagepos(1:N,1), stagepos(1:N,2)]...
%   +[Centerline(1:N,ptAlongWorm,XX), -Centerline(1:N,ptAlongWorm,YY)];

  pos=videoDecimationConversion.*[-stagepos(1:N,1), stagepos(1:N,2)]...
      +[com(1:N,XX), -com(1:N,YY)];

M=1;
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

%quiver(intpos(1:M./2:end-M,1),intpos(1:M./2:end-M,2),diffpos(1:M./2:end,1),diffpos(1:M/2:end,2),'b');
