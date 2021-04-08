videoDecimationConversion=2;

XX=1;
YY=2;
N=iend-istart+1;
Centerline = zeros(N,100,2);
time_speed = zeros(iend-istart+1,1);
for i = 1:N
    j = istart+i-1;
    centerline=reshape(mcd(j).SegmentedCenterline,2,[]);
    Centerline(i,:,1)=centerline(1,:);%x axis
    Centerline(i,:,2)=centerline(2,:);%y axis
    time_speed(i) = mcd(j).TimeElapsed;
end


%% calculte worm postion
%CalculationMethodSelection = 1; %2 means using centroid, 1 means using multiple body points (segmented).
switch CalculationMethodSelection
    case 1 % using multiple points
        segments = 20; %num of body points
        start_counting_segment = 1; % n th point to start calculation
        headpoint = 3; % position intra-segment 
        com = zeros(N,2,segments-start_counting_segment+1);
        for segment = start_counting_segment:segments
        com(:,:,segment-segment+1) = [Centerline(:,headpoint+(segment-1)*100/segments,XX),...
            Centerline(:,headpoint+(segment-1)*100/segments,YY)];
        end
        

    case 2 % using centroid
        segments = 1;
        start_counting_segment = 1;
        com=[mean(Centerline(:,:,XX),2), mean(Centerline(:,:,YY),2)];        
end