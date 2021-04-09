function [rmsp,es,ar,tpred,ts]=predict_nn(xx,y,Tmax,npred,xform,nb,api)
% [rmsp,es,ar,tpred,ts]=predict_nn(xx,tr,Tmax,npred,v,nb,api)
% Nearest neighbor predictor along with estimation of Tpred predict_nn
% error.
% 
% INPUT:
%   xx:    T x m input embedded data (can be any embedding space).
%   y:     Observation time series corresponding to the embeddings
%   Tmax:  Maximum time upto which prediction is to be done (should be long
%          enough so that the error curve saturates). 
%   npred: Number of test indices for prediction.
%   xform: Transformation matrix to project the embedded time-series back
%          to observation space
%   nb:    Number of nearest neighbors
%   api:   range of prediction indices
%        
% OUTPUT:
%   rmsp:   The error curve E(tau) with bootstrapped 95% CI estimates (1 is mean,
%           2 is lower bound, and 3 is upper bound). This scheme is for all
%           estimates.
%   es:     Estimate of the saturation value of the Error curve (loosely
%           speaking the "radius of the attractor", usually related to the
%           variance of the data.
%   ar:     Estimate of the area bounded between es and E(tau).
%   tpred:  Tpred estimate (ar/es).
%   ts:     Estimate of the actual time at which the curve saturates
%
    tdim=size(y,2);
    parfor i=1:npred
        %Get indices of nb (transverse) nearest neighbors
        pred_ind=api(i);
        idx=get_close_recs(xx,pred_ind,nb,Tmax);        
        %average the nb trajectories from idx:idx+Tmax to get the
        %prediction in state space and then project them back to the
        %observation space to calculate the prediction error 
        if ~isempty(idx)                
            nbmat=repmat(idx,1,Tmax+1)+repmat(0:Tmax,length(idx),1);
            yp=squeeze(mean(reshape(xx(nbmat,:),length(idx),Tmax+1,size(xx,2)),1));
            if isrow(yp)
                yp=yp';
            end
            %%yp is the prediction in state space estimated by averaging the neighboring trajectories. Project it back by
            %%transformation v' to get to delay space and then pick the first d columns
            yp=yp*xform';yp=yp(:,1:tdim);   
            yo=y(pred_ind:pred_ind+Tmax,:); %original
            merrp(:,:,i)=((yp-yo).^2); %store the squared error
        else
            continue; 
        end
    end
    %remove points with no neighbors
    zind =squeeze(merrp(1,1,:)==0);   
    merrp(:,:,zind)=[];   
    
    %get bootstrapped estimates (niter bootstraps) of the error curve and related quantities
    %such as the area, es, Tpred and ts.
    niter=200; 
   [rmsp,es,ar,tpred,ts]=process_error_curves(merrp,niter);    
end


function [rmsp,es,ar,tpred,ts]=process_error_curves(err,niter)
    probrmsp=[];
    parfor i=1:niter
        %get a random sample of error curves and average them
        tprmsp(:,i)=sum(sqrt(mean(err(:,:,randi(size(err,3),1,size(err,3))),3)),2);        
        nrms=tprmsp(:,i); 
        
        %get initial estimate of es and mints
        %strictly speaking the following isn't necesarry. It is possible to
        %start with any guess of es and corresponding. The following is a
        %way to get an estimate es which is roughly correct.                
        xl=linspace(min(nrms),max(nrms)+max(std(nrms)),500);   
        [px,xl]=ksdensity(nrms,xl);
        [~,lcs] = findpeaks(px/max(px),xl,'MinPeakHeight',0.5);%[~,pks]=max(pks);
        ines=max(lcs);
        mints=find(nrms>ines,1);

        %estimate of mints by convergence of iterated line fitting
        %the estimate converges in 3-4 iterations. I have set the number of
        %iterations to 10 to be on the safe side.
        ei=(cumsum(nrms));
        for k=1:10            
            tmpfit=abs(polyfit(mints(k):length(ei),ei(mints(k):end)',1));
            mints=[mints;find(nrms>tmpfit(1),1)];            
        end
        
        %check if mints has converged, if not then store those values for
        %debugging
        tmmts=diff(mints);
        cnv(i)=tmmts(end);
        if cnv(i)~=0
            probrmsp=[probrmsp;nrms];        
        end

        %final estimate of paramters
        tmpfit=(polyfit(mints(end):length(ei),ei(mints(end):end)',1));
        tmpes(i)=tmpfit(1);
        tmpar(i)=tmpfit(2);       
        tmpTpred(i)=-tmpfit(2)/(tmpfit(1));
        tmpmints(i)=mints(end);
    end
    %get final estimates and error bars
    rmsp=[sum(sqrt(mean(err,3)),2) prctile(tprmsp,[2.5 97.5],2)];
    es=[mean(tmpes) prctile(tmpes,[2.5 97.5])];
    ar=[mean(tmpar) prctile(tmpar,[2.5 97.5])];
    tpred=[mean(tmpTpred) prctile(tmpTpred,[2.5 97.5])];
    ts=[mean(tmpmints) prctile(tmpmints,[2.5 97.5])];
end