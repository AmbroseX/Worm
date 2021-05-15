function v=wormSpeed(position)
%position是numfram*3数据，三列分别是时间，X，Y坐标 （mm）
vnum=length(position)-1;
v=zeros(vnum,3);
for i=1:vnum
    v(i,1)=position(i,1);
    v(i,2)=(position(i+1,2)-position(i,2))/(position(i+1,1)-position(i,1));  %x速度
    v(i,3)=(position(i+1,3)-position(i,3))/(position(i+1,1)-position(i,1));  %y速度  
end
end