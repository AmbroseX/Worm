function [ flag ] = comparison(x,x1,x2,x3,x4,x5,x6)
%UNTITLED compare slopes 
 if (abs(x-x1)>0.05) && (abs(x-x1)>0.05) && (abs(x-x2)>0.05) && (abs(x-x3)>0.05) && (abs(x-x4)>0.05) && (abs(x-x5)>0.05) && (abs(x-x6)>0.05)
      flag=1;
  else flag=0;

end

