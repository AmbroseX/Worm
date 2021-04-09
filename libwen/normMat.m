function [Theta,normTheta]=normMat(Theta)
    for k=1:size(Theta,2)
        normTheta(k) = norm(Theta(:,k));
        Theta(:,k) = Theta(:,k)/normTheta(k);
    end      
end