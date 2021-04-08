% Given exported analysis data from the Previewer program
% Plot the worm's estimated trajectory as inferred from the stage velocity
% feedback loop.
%
%


load('D:\Analysis\AlkemaCollaboration\flip13\20100916_1653_flp13ChR2Halo6_HUDS_12611-13089_with_stage_vel.mat');

M=10; %smoothing value;

%Transformation between stage velocity space and camera space
T=[-1 0; 0 1;]

%Stage velocity Data
sv=handles.StageVelocity_data*T;

%Raw integrated position
intpos_raw=cumsum(sv);

%Smoothed integrated position
sv_smooth=lowpass1d(sv,M);
intpos=cumsum(sv_smooth);

laserON=T2-T1;
laserOFF=T3-T1;


figure;
hold on;
lw=3;
quiver(intpos(1:M./2:end,1),intpos(1:M./2:end,2),sv(1:M./2:end,1),sv(1:M/2:end,2),'b');
plot(intpos(1:laserON,1),intpos(1:laserON,2),'r','Linewidth',lw);
plot(intpos(laserON:laserOFF,1),intpos(laserON:laserOFF,2),'g','Linewidth',lw);
plot(intpos(laserOFF:end,1),intpos(laserOFF:end,2),'r','Linewidth',lw);
