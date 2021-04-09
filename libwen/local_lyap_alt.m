function [local_lam,ls,nm] = local_lyap_alt(Jmat)
    L=size(Jmat,3);
    At=permute(conj(Jmat),[2 1 3]);
    Aa=cat(3,At,flip(Jmat,3));clear At
    [Q,R,rr]=prod_eig(Aa,2*L);
    M=cat(3,Q(:,:,end),R);
    j=1;
    %lr=log(rr(:,2:end));
    %ls=cumsum(lr(:,k0:end),2)./repmat(1:length(lr(:,k0:end)),size(Jmat,1),1);
    lr(:,j)=sort(infmean(log(rr(:,2:end)),2),'descend');%sort(mean(log(rr(:,2:end)),2),'descend');
    nm(j)=norm(abs(Q(:,:,end))-eye(size(Jmat,1)));
    Qe(:,:,j)=Q(:,:,end);
    while(nm(j)>0.0001 && j<100)
        j=j+1;
        [Q,R,rr]=prod_eig(M,2*L+1);  
        M=cat(3,Q(:,:,end),R(:,:,2:end));
        lr(:,j)=sort(infmean(log(rr(:,2:end)),2),'descend');%sort(mean(log(rr(:,2:end)),2),'descend');
        nm(j)=norm(abs(Q(:,:,end))-eye(size(Jmat,1)));    
        Qe(:,:,j)=Q(:,:,end);
    end
    %local_lam=lr(:,amin(nm));
    if j<100
        local_lam=lr(:,end)';
        Qe=Qe(:,:,end);
        nm=nm(end);
    else
        local_lam=median(lr');
        Qe=Qe(:,:,amin(nm));
        nm=nm(amin(nm));
    end    
  
    lr=log(rr(:,2:end));
    lr(isinf(lr))=nan;
    ls=cumsum(lr(:,1:end),2,'omitnan')./repmat(1:length(lr(:,1:end)),size(Jmat,1),1);
    ls=(downsample(ls',2));
end