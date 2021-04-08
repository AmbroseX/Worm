function r=estimate_turning_rate(y,t)
[~,locs] = findpeaks(-y,'MINPEAKHEIGHT',-75,'MINPEAKDISTANCE',10);
sigma_w=2;

N=length(t);

r=zeros(N,1);

for i=1:length(t)
    r(i)=0;
    for j=1:length(locs)
        r(i)=r(i)+1/sqrt(2*pi)/sigma_w*exp(-(t(i)-t(locs(j)))^2/2/sigma_w^2); 
    end
end



        
