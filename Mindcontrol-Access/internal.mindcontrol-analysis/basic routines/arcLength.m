function s=arcLength(curvx, curvy)
% Function to find the arc length
%
% Andrew Leifer
% leifer@fas.harvard.edu
% 1 November 2010

dx=curvx(1:end-1)-curvx(2:end);
dy=curvy(1:end-1)-curvy(2:end);

s=sum([dx.^2+dy.^2].^0.5);

end