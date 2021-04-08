body_cur=zeros(numframes,99);
body_amp=zeros(1,99);
for bodyseg=2:99
    
    origin=bodyseg;
    radius=1;
    for j=1:numframes
        body_cur(j,bodyseg)=mean(curvdatafiltered(j,origin-radius:origin+radius));
    end
end
for bodyseg=1:99
    body_amp(1,bodyseg)=std(body_cur(:,bodyseg));
end
disp(body_amp);