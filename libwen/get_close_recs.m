function [idx,disi]=get_close_recs(xx,pind,nb,Tmax,eps)
% [idx,disi]=get_close_recs(xx,pind,nb,Tmax)
% Find nb close (transverse) recurrences to the point pind in the state space xx.
% 
% INPUT:
%   xx:      T x m input embedded data (can be any embedding space).
%   pind:    Query point whose neighbors have to be determined.
%   nb:      Number of nearest neighbors asked
%   Tmax:    (optional) used when length(xx)-Tmax is an upper limit to the indices.
%   eps: (optional) only return neighbors with distance less than eps
%        
% OUTPUT:
%   idx:     Time indices of the recurrences.
%   disi:    Distances corresponding to the time indices.


%%Estimate distances from point indexed by pind to every other point
    dist=sqrt(sum(bsxfun(@minus,xx,xx(pind,:)).^2,2));dist(dist==0)=nan; %get distance     
%%Find close recurrences by identifying local minima of the distance function.
%%Flip so that minima are maxima; Normalize between 0 and 1 and set a 5% prominence threshold to get rid of
%%small peaks.  
    ndist=norm2one(-dist);
    [pks,lcs]=findpeaks(norm2one(-dist),'MinPeakProminence',0.05); 
%%Findpeaks doesn'take the endpoints into account, include it explicitly.
    pks=[dist(lcs);dist(1)];lcs=[lcs;1];  
%%If Tmax has been set then remove all indices that get out of bounds when
%%propagated
    if exist('Tmax','var')
        pks((lcs+Tmax>length(xx)))=[];
        lcs((lcs+Tmax>length(xx)))=[];
    end
    
   
    
    if isempty(pks)        
        idx=[];
        disi=[];
    else        
%%sort indices in the order of nearest neighbor first
        idx=asort(lcs,pks); 
        disi=dist(idx);     
        
        if exist('eps','var')
            idx=idx((disi<eps));
            disi=disi((disi<eps));
        end
        
        if length(idx)>nb           
            idx=idx(1:nb);
            disi=disi(1:nb);%*mean(dist,'omitnan');
        end        
    end  
end