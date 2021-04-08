function [counts,binned_t]=count_turns(y,t)
[~,locs] = findpeaks(-y,'MINPEAKHEIGHT',-75,'MINPEAKDISTANCE',40);

bin_size=1; 

binned_t=0:bin_size:t(end);

t_peaks=t(locs);

counts=histc(t_peaks,binned_t);

