function [ampmatb, ampmatd]=extractampDorB(amp)
ampmatb=zeros(length(amp),100);
ampmatd=zeros(length(amp),100);
for i=1:length(amp)
    ampmatb(i,:)=amp{i,1}(1,:);
    ampmatd(i,:)=amp{i,1}(2,:);
end
end