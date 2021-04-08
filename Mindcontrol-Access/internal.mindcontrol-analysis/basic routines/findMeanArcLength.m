function s=findMeanArcLength(centerline,sigma)
% Finds the mean arc lenght of a 3D matrix, centerline, where
% the first dimension is the instance to be averaged over
% the second dimension is the pt along the centerline
% and third dimension is the x or y coordinate respectively
%
% For example: 
% size(centerline)
% ans =
%     318 100 2
% refers to the case where there are 318 frames, with 100 points along the
% centerline
% 
%
% Sigma specifies the rise time of a lowpass filter 
%
% by Andrew Leifer
% 1 November 2010
% leifer@fas.harvard.edu


for k=1:size(centerline,1)
    c(k)=arcLength(lowpass1d(centerline(k,:,1),sigma),lowpass1d(centerline(k,:,2),sigma));
end
s=mean(c);
end