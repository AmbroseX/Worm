function center=centeroidPosion(A)
time=A.TimeElapsed;
X=A.Centerline(:,:,1);
Y=A.Centerline(:,:,2);

Xc=mean(X,2);
Yc=mean(Y,2);
numfram=length(Xc);
center=zeros(numfram,3);
center(:,1)=time;
center(:,2)=Xc;
center(:,3)=Yc;

end