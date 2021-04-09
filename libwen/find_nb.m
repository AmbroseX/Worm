function Nb=find_nb(xx,ntind,nbrange)
% Nb=find_nb(xx,ntind,nbrange)
% Find number of nearest neighbors by minimizing the 1 step prediction
% error.
% 
% INPUT:
%   xx:      T x m input embedded data (can be any embedding space).
%   ntind:   Number of test indices for prediction.
%   nbrange: Range of neighborhood indices to search in.
%        
% OUTPUT:
%   Nb:     Number of nearest neighbors.

    len=size(xx,1);
    if ~exist('ntind','var')
        ntind=floor(len/10);
    end
    if ~exist('nbrange','var')
        nbrange=50;
    end
    pind=randi(len-1,1,ntind); %get ntind random test points
    dd=ones(ntind,nbrange)*nan; %initialize error matrix
    parfor ppi=1:length(pind)
        %calculate distances and set self distance to nan       
        dist=sqrt(sum(bsxfun(@minus,xx,xx(pind(ppi),:)).^2,2)); dist(dist==0)=nan; 
        %Find close recurrences by identifying local minima of the distance function.
        %Flip so that minima are maxima; Normalize between 0 and 1 and set a 5% prominence threshold to get rid of
        %small peaks.  
        dist=norm2one(-dist);
        [pks,lcs]=findpeaks(dist,'MinPeakProminence',0.05); 
        %findpeaks doesn'take the endpoints into account, include it explicitly.
        pks=[pks;dist(1)];lcs=[lcs;1]; 
        %sort indices in the order of nearest neighbor first
        lcs=asort(lcs,pks,'descend'); 
       
        %Find the 1 step prediction error for each set of NN values
        for i=1:nbrange 
            if i>length(lcs)
                break;
            else              
                dd(ppi,i)=sum(abs(mean(xx(lcs(1:i)+1,:))-xx(pind(ppi)+1,:))); %error in predicting the (i+1)th point from its nearest neighbors
            end
        end
    end 
    
    %get bootstrapped mean error as a function of nearest neighbors
    opt = statset('UseParallel',true);
    stats=bootstrp(1000,@mean,bsxfun(@rdivide,dd',min(dd'))','Options', opt);
    lec=[mean(bsxfun(@rdivide,dd',min(dd'))')' prctile(stats,[2.5 97.5])'];
    
    %get bootstrapped error variance as a function of nearest neighbors
    opt = statset('UseParallel',true);
    stats=bootstrp(1000,@std,bsxfun(@rdivide,dd',min(dd'))','Options', opt);
    lecst=[std(bsxfun(@rdivide,dd',min(dd'))')' prctile(stats,[2.5 97.5])'];
    
    %We want nearest neighbor values with low variance and low bias
    biasvar=sqrt(lecst.^2+lec.^2);
    Nb=max(amin(biasvar));   
end