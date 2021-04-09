% Jmat=cat(3,eye(size(Jmat,1)),Jmat);
% Jmat = Jmat(:,:,1:end-1);
% R(:,:,1)=eye(size(Jmat,1),size(Jmat,1));
% Q(:,:,1)=eye(size(Jmat,1),size(Jmat,1));
% rr(:,1)=diag(R(:,:,1));
% 
% for i=2:length(Jmat);
%     %[Q(:,:,i), R(:,:,i)]=gsorth(Jmat(:,:,i)*Q(:,:,i-1));
%     %rr(:,i)=diag(R(:,:,i));
%     [Q(:,:,i), R(:,:,i)]=qr(Jmat(:,:,i)*Q(:,:,i-1));
%     rr(:,i-1)=abs(diag(R(:,:,i)));
%     disp(i);
% end
% lr=log(rr);
% lr = lr(:,~isnan(lr(1,:)));
% ls=cumsum(lr(:,200:end),2)./repmat(1:length(lr(:,200:end)),size(Jmat,1),1);
% 

function [Q,R,rr]=prod_eig(seqMat,L)
    seqMat=cat(3,eye(size(seqMat,1)),seqMat);
    R(:,:,1)=eye(size(seqMat,1),size(seqMat,1));
    Q(:,:,1)=eye(size(seqMat,1),size(seqMat,1));
    rr(:,1)=diag(R(:,:,1));
    for i=2:L+1
%         [Q(:,:,i), R(:,:,i)]=gsorth(seqMat(:,:,i)*Q(:,:,i-1));
%         rr(:,i)=diag(R(:,:,i));
        [Q(:,:,i), R(:,:,i)]=qr(seqMat(:,:,i)*Q(:,:,i-1));
        rr(:,i-1)=abs(diag(R(:,:,i)));        
    end
    Q=Q(:,:,2:end);
    R=R(:,:,2:end);  
    %rr=rr(:,1:end);
end