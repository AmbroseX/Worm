function [merr,Jac_data]=learn_jac(xx,T1,tmax,rmax,eps)
% [merr,Jac_data]=learn_jac(xx,T1,tmax,rmax,eps)
% Detect periodic orbits by finding close recurrences.% 
% 
% INPUT:
%   xx:   T x m input embedded data (can be any embedding space).
%   T1:   Estimate T1-step Jacobian (default is 1, but can be changed.         
%   tmax: The maximum time index to search for a periodic orbit. tmax/pmin
%         is the maximum period sequence which will be output.
%   rmax: Number of recurrences used to calculate the function
%         epsilon(r,tau)'. Default is half the number of time-points.
%   eps : Epsilon used for the weighting scheme in Jacobian estimation.
%         By default eps is calculated from the recurrence function.
%
% OUTPUT:
%   merr: prediction error for the local Jacobians
%   Jac_data: Contains the Jacobians and lyapunov spectrum estimates.
  %% Set defaults
    m=size(xx,2);
    if ~exist('T1','var')
        T1=1;
    end
    
    if ~exist('tmax','var')
        tmax=750;
    end
    if ~exist('rmax','var')
        rmax=floor(length(xx)/2);
    end
  
    %% prepare the input and output matrices. 
    %Input is the current state, output is the state after T2 steps
    XX=xx(1:end-T1,:);
    Xn=xx(T1+1:end,:);
    
    %% Find the minimum scale eps. 
    if ~exist('eps','var')
        parfor i=1:tmax
            tp=sqrt(sum((xx(1:end-i,:)-xx(i+1:end,:)).^2,2));    
            tp=sort(tp);
            for j=1:rmax
                epsrtau(j,i)=tp(j);
            end
        end
        mn=mean(epsrtau);
        warning('off','signal:findpeaks:largeMinPeakHeight')
        [~,pers,~,~]=myfindpeaks(norm2one(-mn),'Annotate','extents','WidthReference','halfheight','minpeakprominence',0.01);    
        epstau=sgolayfilt(epsrtau(m,:),1,7);  
        eps=epstau(pers(1));
    end
    
%%    
    lib=1:length(XX);   
    parfor ipred=1:length(XX)

        %now fit on all space
        libs=lib;libs(ipred)=[];
        dist=pdist2(XX(libs,:),XX(ipred,:));

        wts=exp(-dist/eps);

        %give zero weight to successive points, so they don't bias the
        %least square fit.
        if ipred<=15
            wts(1:15)=0;
        elseif (length(dist)-ipred)<=7-1
            wts(ipred-7:end)=0;
        else
            wts(ipred-7:ipred+7)=0;
        end

       
        % the following removes entries of small weights. It helps speed up
        % the estimation process, but could lead to inaccurate estimates.
        %libs((wts/max(wts))<1e-3)=[];  
        %wts((wts/max(wts))<1e-3)=[];
        
        z=[ones(length(libs),1) XX(libs,:)];z=bsxfun(@times,z,wts);
        zn = bsxfun(@times,Xn(libs,:),wts);
        A=z\zn;

        AA(:,:,ipred)=A;

        %2 steps error in prediction used only to judge the quality of the
        %jacobian estimates.
%         T2=2;
%         if (ipred+T2)>length(XX)
%             T2=length(XX)-ipred;
%         end        
%         xi=XX(ipred,:);
%         for ni=1:T2+1
%             xi=[xi;[1 XX(ipred,:)]*A];
%         end
%         err(ipred,:)=(Xn(ipred+T2,:)-xi(end,:)).^2;

            
        xi=[1 XX(ipred,:)]*A;
        err(ipred,:)=(Xn(ipred,:)-xi).^2;

        %%This is another method of estimating Jacobians which is common.
        %%Instead of fitting sequential dynamics, we first approximate the
        %%tangent space by the difference vectors and then fit a linear
        %%system to the difference vectors. In principle, Amat and Jmat
        %%should be the same, but in practice there are differences.
        z = bsxfun(@minus,XX(libs,:),XX(ipred,:));
        z = bsxfun(@times,z,wts);
        zn = bsxfun(@minus,Xn(libs,:),Xn(ipred,:));
        zn = bsxfun(@times,zn,wts);    
        Jmat(:,:,ipred)=(pinv(z)*zn)';
    end            
    opt = statset('UseParallel',true);
    stats=bootstrp(1000,@(x) sum(sqrt(mean(x))),err,'Options', opt);
    merr=[sum(sqrt(mean(err))) prctile(stats,[2.5 97.5])];


    %%store jac data   
    Jac_data.eps=eps;
    %Jac_data.per=per;
    Jac_data.Amat=AA;
    Jac_data.Jmat=Jmat;
    
    %estimate lyapunov exponents by recursive QR
    [local_lam,ls]=local_lyap_alt(AA(2:end,:,:));
     Jac_data.lexp=local_lam;
     Jac_data.lexp_avg=ls;
end