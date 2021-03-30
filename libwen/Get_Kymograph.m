%This function Plot Kymograph of Curvature by the Sequence of Modes
%and Get Backbone for all Modes:
%
%V is the embedding result of svd 
%eigenvectors is the PCA results of curvature
%Typically: K_mode=5, embedded_dimension=6
%
%
%28 October 2019
%Jiahao Tang
%tangjh0513@163.com




function  [Kymo_curvdata,Kymo_theta] = Get_Kymograph(Vkymo,RatioKymo2Curvdata,eigenvectors,K_delay,K_mode,embedded_dimension)
    %Reshape the embedding vectors V to V_reshaped:
    V_reshaped = reshape(Vkymo,[K_mode,K_delay,embedded_dimension]);
    %Define new variables:
    Kymo_curvdata = zeros(size(eigenvectors,1)-1,size(V_reshaped,2),size(V_reshaped,3));
    Kymo_theta = zeros(size(eigenvectors,1),size(V_reshaped,2),size(V_reshaped,3));
    Kymo_backbone_x = zeros(size(eigenvectors,1),size(V_reshaped,2),size(V_reshaped,3));
    Kymo_backbone_y = zeros(size(eigenvectors,1),size(V_reshaped,2),size(V_reshaped,3));

    %RatioKymo2Curvdata is the scalar factor set to get curvature in mm^-1, 
    %it is not accurate and may need adjustment

    for m = 1:embedded_dimension
        %Get curvature
        Kymo_theta(:,:,m) = eigenvectors(:,1:K_mode)*V_reshaped(:,:,m);
        Kymo_curvdata(:,:,m) = (RatioKymo2Curvdata).*diff(Kymo_theta(:,:,m),1,1);
        %Plot kymograph of curvature
        %Get backbone
        %Warning: the backbone is not completely accurate for now
        %It depends greatly on RatioKymo2Curvdata
        Kymo_theta(:,:,m) = Kymo_theta(:,:,m) - repmat(mean(Kymo_theta(:,:,m),1),100,1);    
        Kymo_backbone_x(:,:,m) = cumsum(cos((pi/180).*Kymo_theta(:,:,m)));
        Kymo_backbone_y(:,:,m) = cumsum(sin((pi/180).*Kymo_theta(:,:,m)));   
    end
    
    %lim = [min(min(min(Kymo_curvdata))), max(max(max(Kymo_curvdata)))];
    figure;
    for m = 1:embedded_dimension
    subplot(ceil(embedded_dimension/2),2,m);imagesc(Kymo_curvdata(:,:,m)',[-5,5]);title(['Mode' num2str(m)]);colorbar;colormap jet;
    end


%If you want to see the backbone,run this loop:
%Warning:this part is still being tested, this backbone is not accurate, and depends greatly on RatioKymo2Curvdata
%{


p=1;    % you may change p, p= 1:embedded_dimension,
for q = 1:K_delay
    %Plot backbone one by one for given embedded_dimension
    figure;plot(Kymo_backbone_x(:,q,p),Kymo_backbone_y(:,q,p))    
end   


%}



end
