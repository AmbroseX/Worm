%By Jiahao Tang September 30, 2019
%For questions, email tangjh0513@163.com
%this function allow you to create a 3D trajectory(x,y,z) with speed s
%s=0,1,2,....
%You can copy this to the end of your code, 
%or add the animation3.m to your folder
function animation3(x,y,z,s)
    
    h = animatedline('Color',[0 0.4470 0.7410],'LineWidth',1);
    axis([min(x),max(x),min(y),max(y),min(z),max(z)])
    
    for k = 1:s+1:length(x)-s
        xvec = x(k:k+s);
        yvec = y(k:k+s);
        zvec = z(k:k+s);
        addpoints(h,xvec,yvec,zvec)
        drawnow
    end  
end