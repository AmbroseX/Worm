tic
%% load data
clear,clc
% N2 silence testDLP
filepath='N2';

%for Rongkang desktop-3070  & Laptap
workpath=fullfile('G:','Data','WenLab','Worm_Embed');
%For the 2080Ti
% workpath=fullfile('/','home','wenlab','xrk','Worm_Embed');

addpath(genpath(fullfile(workpath,'libwen')));
addpath(genpath(fullfile(workpath,'data',filepath)));

disp('Staring load data...')

load('20210429_2106_w3.mat') % 1*12 cell ,33600*5 double


%calulate the relative speed.
position=wormrelativePosion(wormdata);
speed=wormSpeed(position);