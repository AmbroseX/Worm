function [ difflag ] = judge( slope,sloperef,setpoint,j1,jmin2,j,i )
%Judge whether the curvature is apparently different from the referee
%Using slope 
totalbias=[0,0,0,0,0];
        for k=1:5
            totalbias(i) = totalbias[i]+(sloperef(j+j1-jmin2(i)+k-1,i)-slope(j+k-1,i))^2
        end
 if difflag && (sloperef(j+j1-jmin2(i),i)-slope(j,i))^2 > setpoint && totalbias(i)/5 > setpoint
            difflag = 1;
 else difflag = 0;
     
 end
