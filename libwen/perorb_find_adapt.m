function [pertab,pers,epstau]=perorb_find_adapt(xx,tmax,rmax,rupo)
% [finrec,pers,in_d]=perorb_find_adapt(xx,tmax,rmax,rupo)
% Detect periodic orbits by finding close recurrences.% 
% 
% INPUT:
%   xx:   T x m input embedded data (can be any embedding space).
%   tmax: The maximum time index to search for a periodic orbit. tmax/pmin
%         is the maximum period sequence which will be output.
%   rmax: Number of recurrences used to calculate the function
%         epsilon(r,tau)'. Default is half the number of time-points.
%   rupo: Number of recurrences used to find the minimum epsilon. Default
%         is rupo=m, the number of embedding dimensions.
%
% OUTPUT:
%   pertab: Table of periodic orbits.  [time_index, distance_of_approach, number_of_continuous_elements, length, angle, integer_period]
%   pers  :   Period lengths estimated by local minima of epstau.
%   epstau: Filtered epsilon(r,tau) for r=rupo. This sets the epsilon to
%           find recurrences

%% Set defaults
    m=size(xx,2);
    if ~exist('tmax','var')
        tmax=750;
    end
    if ~exist('rmax','var')
        rmax=floor(length(xx)/2);
    end
    
     if ~exist('rupo','var')
        rupo=m;
    end
   
%% Calculate the epsilon-tau function. Get the periods, their widths and epsilon.
    parfor i=1:tmax
        tp=sqrt(sum((xx(1:end-i,:)-xx(i+1:end,:)).^2,2));    
        tp=sort(tp);
        for j=1:rmax
            epsrtau(j,i)=tp(j);
        end
    end
    mn=mean(epsrtau);
    warning('off','signal:findpeaks:largeMinPeakHeight')
    [~,pers,~,wx]=myfindpeaks(norm2one(-mn),'Annotate','extents','WidthReference','halfheight','minpeakprominence',0.01);    
    epstau=sgolayfilt(epsrtau(rupo,:),1,7);  
    eps=max(epstau(pers));
%% Get the recurrences and their statistics
    parfor i=1:size(xx,1)
        dist=sqrt(sum(bsxfun(@minus,xx,xx(i,:)).^2,2));dist(dist==0)=nan; %find distances from xx(i) to every point
        warning('off','signal:findpeaks:largeMinPeakHeight')
        [pks,lcs]=findpeaks(-dist,'MinPeakHeight',-eps); % find all peaks with distance less than eps
        %only consider future neighbors
        pks=pks(lcs>i);
        lcs=lcs(lcs>i); 
        
        %extract statistics
        if ~isempty(lcs)
            rec(i)=lcs(1)-i; %time for first time state recurs (gets within eps of x(i)).
            recdis(i)=-pks(1); %distance for first recurrence                        
        else
             rec(i)=nan;
             recdis(i)=nan;           
        end
    end

%% Process the recurrence statistics above and store them for periods upto tmax
    
    jind=1:tmax; % min to max period
    bigst=[];
    idx=1;
    for k=1:length(jind)
        ti=find(rec==jind(k)); %index of the recurrence of periods of length jind(k)
        if isempty(ti),continue;end %if there's no recurrence
        rdi=recdis(ti); %distance of the recurrence   

        bi=[1 find(diff(ti)>1)+1]'; %find successive indices
        sti=[bi [bi(2:end)-1;length(ti)]]; %get starting and ending indices of the succesive indices in previous step
        for ii=1:size(sti,1)
            [~,am]=min(rdi(sti(ii,1):sti(ii,2)));am=am+sti(ii,1)-1; %get the index of the closest recurrence            
            %get the indices of all points on the periodic orbit
            idx=ti(am):ti(am)+jind(k)-1;             
            %get the angle between first and last point on the periodic
            %orbit. Again used to weed out spurious orbits as the angle
            %should be the close. A better solution is required to detect
            %and remove spurious recurrences.
            if length(idx)>3
                p1=(xx(idx(1)+1,:)-xx(idx(1),:));
                p2=(xx(idx(end)+1,:)-xx(idx(end),:));
                ng=dot(p1/norm(p1),p2/norm(p2));
            else
                ng=nan;
            end
            %store all the statistics in the following order
            %[time_index, distance_of_approach,
            %number_of_continuous_elements, period, angle
            sten(ii,:)=[ti(am) rdi(am) sti(ii,2)-sti(ii,1)+1 jind(k) ng];
        end
        bigst=[bigst;sten];clear sten;
    end
    
    %keep only those orbits where the first and last point are nearly
    %parallel. Needs a better solution.
    sel=bigst(bigst(:,5)>0.9,:);    
    
%% Calculate the integer period using the widths of the peaks wx estimated above    
    pertab=[];
    for pki=1:size(wx,1)
        sbig=sel(find(sel(:,4)>floor(wx(pki,1)) & sel(:,4)<floor(wx(pki,2))),:); %select all orbits within a peak
        if isempty(sbig),continue,end
        sbig=asort(sbig,sbig(:,1)); %sort by time index

        %remove duplicate orbits, i.e. orbits that are close to each other
        %in time but have the same integer period.
        flag=0;
        while(1)
            flag=0;
            ti=1;
            while (flag==0)
               
                tti=find(sbig(:,1)<=sbig(ti,1)+pers(pki));
                if tti(end)==ti
                    flag=1;
                else
                    ti=tti(end);
                end
            end
            mi=amin(sbig(tti,2));
            pertab=[pertab;[sbig(mi,:) pki epstau(pers(pki))]];
            sbig(tti,:)=[];
            if isempty(sbig),break;end
        end
    end
   pertab(pertab(:,2)>pertab(:,7),:)=[];
   pertab(:,7)=[];    
%   finrec=getTimeScale(xx,finrec,eps);
   warning('on','signal:findpeaks:largeMinPeakHeight')
end

%%
function finrec=getTimeScale(xx,finrec,eps)

    for j=1:size(finrec,1)
        idx=finrec(j,1):finrec(j,1)+finrec(j,4); %get index of the UPO points
        alni=[];
        parfor i=1:length(idx)
        [ni,dis]=get_close_recs(xx,idx(i),500,2500,eps);
        alni=[alni;ni];
        end
        alni=unique(alni);
        %propogate them in time and count how many remain within epsilon
        nn=[];
        parfor i=1:2500
        dd=pdist2(xx(idx,:),xx(alni+i-1,:));
        nn(i)=length(find(min(dd)<eps));
        end
        tt=((1./(1:length(nn))).*abs(log(nn/nn(1))));
        wt=(nn/nn(1)).^(finrec(j,4)./(1:length(nn)));
        tau(j)=wt(end);  
        disp(j);
    end  
    
    finrec=[finrec tau'];
end