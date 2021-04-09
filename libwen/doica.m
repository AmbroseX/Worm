function [sicasig, A, W, white, dewhite, whitesig]=doica(input,numcomps,numit)
% [sicasig, A, W, white, dewhite, whitesig]=doica(input,numcomps,numit)
% Wrapper function for fastica
% 
% INPUT:
%   input:   Kd x T observation space.
%   numcomps: number of independent components required (this should equal
%   the embedding dimensiom estimated using Tpred calculation
%   numit: Number of iterations.

%
% OUTPUT:
%   icasig: Phase space trajectories
%   A  :   Mixing Matrix (behavioral modes)
%   W  :   Demixing Matrix
%   white: Whitening Matrix
%   dewhite: Dewhitening Matrix
%   whitesig: Whitened Signal

%%
    if nargin<3
       numit=10000;
    end
    k=1;
    while 1
        if k>1, break;end
        k=k+1;
        [icasig, A, W, white, dewhite, whitesig]=fastica(input,'lastEig',numcomps,'approach','symm', 'epsilon',1e-7,'maxNumIterations',numit,'g','pow3','finetune','pow3','stabilization','on','maxFinetune',500);
        if size(icasig,1)==numcomps, break;end
    end
    if size(icasig,1)~=numcomps, sicasig=[];A=[];W=[]; return;end
    
    %sort components by reconstruction error  (Not essential, just a
    %reordering).
    for i=1:numcomps
        inds=setdiff(1:numcomps,i);
        e=A(:,inds)*icasig(inds,:) - input;
        fe(:,i)=sum(e.^2,2)/length(e);
        disp(i)
    end
    [sfe, ind]=sort(sum(fe),'descend');
    sicasig=icasig(ind,:);
    A=A(:,ind);
    W=W(ind,:);
end