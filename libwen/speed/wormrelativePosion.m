function relativePosion=wormrelativePosion(A)
%算质心速度speed随时间变化的数据
micronperunit = 0.05;  %每个单位位移代表的长度 mm
micronperpxl = 2.43;  %每个像素代表的长度 mm
time=A.TimeElapsed;
centerPosion=centeroidPosition(A);   %质心相对相机的时间-位移 X1

StagePosition=A.StagePosition;
snumfram=length(time);

relativePosion=zeros(snumfram,3);

for i=1:snumram
   relativePosion(i,1)=time(i);
   relativePosion(i,2)=micronperpxl*centerPosion(i,2)-micronperunit*StagePosition(i,2);
   relativePosion(i,3)=micronperpxl*centerPosion(i,3)-micronperunit*StagePosition(i,3);
end
end