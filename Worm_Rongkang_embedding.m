tic
%% load data
clear,clc
% N2 silence testDLP
filepath='testDLP';

%for Rongkang desktop-3070  & Laptap
workpath=fullfile('G:','Data','WenLab','Worm_Embed');
%For the 2080Ti
% workpath=fullfile('/','home','wenlab','xrk','Worm_Embed');

addpath(genpath(fullfile(workpath,'libwen')));
addpath(genpath(fullfile(workpath,'data',filepath)));

disp('Staring load data...')

load('20210403_2042_w3.mat') % 1*12 cell ,33600*5 double


wormName = wormdata.wormname;   %to create folder to keep .jpg
savefolder=fullfile(workpath,'prodata',filepath,wormName);

if exist(savefolder)==0
    disp('dir is not exist');
    mkdir(savefolder);
    disp('make dir success');
else
    disp('dir is exist');
end

%% Our data: preparation
angle_data=wormdata.angle_data;
curve_data=wormdata.curv_data;

disp('loadding success ...')

ts=16; %sampling rate
rtheta_s = angle_data(:,2:end); %去掉第一列 得100列维度
mean_X = mean(rtheta_s, 1); %整列做平均
T_sample = size(rtheta_s, 1);
X_shifted = rtheta_s - repmat(mean_X, T_sample, 1); %将rtheta_s减去均值，中心化

K_mode = 5; %
r = 4; %ratio of folding
T = ceil(T_sample / r); %downsampled data points
L = size(rtheta_s, 2); %number of segments 
X_downsampled = zeros(T, L);

for i = 1:L
    X_downsampled(:, i) = decimate(X_shifted(:, i), r); %decimate向下采样时会过滤掉再进行采样，防止混叠
end

%% Our data: get eigenworms by PCA
[eigenvectors, eigenvalues] = eig(X_downsampled' * X_downsampled);
eigenvectors = fliplr(eigenvectors); %如果 A 是一个行向量，则 fliplr(A) 返回一个相同长度的向量，
%其元素的顺序颠倒。如果 A 是一个列向量，则 fliplr(A) 只返回 A，对于多维数组，fliplr
%作用于由第一个和第二个维度构成的平面。特征值由大到小的排列的向量
c = (X_downsampled * eigenvectors);  %投影到eigenvectors张成的空间
y = c(:, 1:K_mode); %the coefficients of the first K_modes

%% optimize K
%params
%Kmax, maximum K used in the search. As a first guess, autocorrelation
%time, zero of autocorrelation function,oscillation period or its double.
%Delrange can be used to sample K coarsely, recommended when starting.
maxdel = 80; delrange = 1:1:maxdel;
%number of random points used for prediction. Larger number of points give
%a better estimate, but takes longer.
npred = 7000; %7000
%tau_max. Harder to estimate, set it to several multiples of the
%autocorrelation function. Check and error curves in rmsp_d and verify that
%the error has saturated. If it never saturates then most likely the
%nonstationarity is too strong. Contact tosifahamed@gmail.com in this case
Tmax = 700;

%generate random test indices for the search
[ybar, ~, idx, ~] = embed_cells(y, maxdel); %embed for maximum delay to calculate the valid range of test indices
%center ybar by removing mean
ybar = bsxfun(@minus, ybar, mean(ybar));
nb = find_nb(ybar); %find nearest neighbors and use the same for the entire search
%getting valid test indices
idx = idx(1:end - Tmax, 1);
testinds = idx(randi(length(idx), npred * 2, 1));

%search for K
for deli = 1:length(delrange)
    [ybar, ybar_orig, ybar_orig_idx] = embed_cells(y, delrange(deli));
    %get test indices that are included in the embedding
    %nidx=embed_cells({oi},delrange(deli));nidx=nidx(:,1);
    [~, tstidx] = ismember(testinds, ybar_orig_idx);
    [rmsp_k{deli}, es_k(deli, :), ar_k(deli, :), tpred_k(deli, :), ts_k(deli, :)] = predict_nn(ybar, ybar_orig, Tmax, npred, 1, nb, tstidx);
    disp(delrange(deli));
end


%plot Tpred(k)
figure
hold on
errorbar(delrange, tpred_k(:, 1), tpred_k(:, 1) - tpred_k(:, 2), tpred_k(:, 3) - tpred_k(:, 1), 'ko', 'markersize', 8);
%errorbar(x,average,'Var_nagetive','Var_positive');
%title('Search for optimal K')
xlabel('K')
ylabel('T_{pred}')
hold off
saveas(gcf, fullfile(savefolder,'Tpred(K).jpg'))

[pks,locs] = findpeaks(tpred_k(:,1));
Kstar=locs(end);  % K越大越好
Kstar=41;
disp(['Best K is:',num2str(Kstar)]);

%% Use SVD to optimize m
%set the value of Kstar as determined from Tpred(K) above
%Maximum embedding dimension corresponds to the maximum number of singular
%vectors to consider.
[ybar, ybar_orig, ybar_orig_idx] = embed_cells(y, Kstar);
[~, tstidx] = ismember(testinds, ybar_orig_idx);
%perform SVD，奇异值分解寻找最佳的m
[U, S, V] = svd(ybar, 'econ');

mmax = 14; %m最大到多少为止
for mi = 1:mmax
    [rmsp_m{mi}, es_m(mi, :), ar_m(mi, :), tpred_m(mi, :), ts_m(mi, :)] = predict_nn(U(:, 1:mi), ybar_orig, Tmax, npred, V(:, 1:mi) * S(1:mi, 1:mi), nb, tstidx);
    disp(mi);
end
% plot Tpred(m)
figure
hold on
errorbar(1:mmax, tpred_m(:, 1), tpred_m(:, 1) - tpred_m(:, 2), tpred_m(:, 3) - tpred_m(:, 1), 'ko', 'markersize', 8);
xlabel('m')
ylabel('T_{pred}')
hold off

saveas(gcf, fullfile(savefolder,'Tpred(m).jpg'));

mStar = 6;  % to set best embedding dimension

%center ybar by removing mean
ybar = bsxfun(@minus, ybar, mean(ybar));
Y = ybar;
%% Get Delay matrix by Wen
Kymoratio_Jiahao = 100;
K_delay = 12;
K_mode_dimension = 5;
embedded_dimension = mStar;
num_Frame = size(y, 1);   %y是前面K_modes个modes的矩阵系数
Y = zeros(num_Frame - K_delay + 1, K_mode_dimension * K_delay); %初始化embedding后的矩阵大小

%embedding the data
j = 1;
for k = 1:K_delay
    Y(:, j:j + K_mode_dimension - 1) = y(k:k + num_Frame - K_delay, :);
    j = j + K_mode_dimension;
end
ybar=Y;

%% Get modes&trajectory by SVD or ICA
%using covairance matrix
%[V,D] = eig(Y'*Y);
%V = fliplr(V);
%c_embedded = Y*V;

Y(isnan(Y)) = 0; %让Y中所有为nan的元素为0
[U, S, V] = svd(Y, 'econ');
XX = (Y * V(:, 1:embedded_dimension))';
len = floor(size(XX, 2) / 2);

figure
hold on
title('SVD phase plot')
plot3(XX(1, end - len:end), XX(2, end - len:end), XX(4, end - len:end));
xlabel('X_t_1')
ylabel('X_t_2')
zlabel('X_f_2')
hold off
saveas(gcf,fullfile(savefolder,'SVD_phase.jpg'))




%% ICA by our method:Centered and whitened
%mu = mean(X',2);
%T = sqrtm(inv(cov(X)));
%Zcw = T * bsxfun(@minus,X',mu);
mu = mean(Y);
[X, whiteningMatrix, dewhiteningMatrix] = whitenRows(bsxfun(@minus, Y, mu)', embedded_dimension);
[~, B] = fastICA(X, embedded_dimension, 'negentropy');
deMix = B' * whiteningMatrix;  %白化矩阵
mix = dewhiteningMatrix * B;
Zica = deMix * Y' + (deMix * mu') * ones(1, size(Y, 1));
[Kymo_curvdata_Jiahao, Kymo_theta_Jiahao] = Get_Kymograph(mix, Kymoratio_Jiahao, eigenvectors, K_delay, K_mode, embedded_dimension);

%icasig = W * mixedsig + (W * mu) * ones(1, NumOfSampl);
%check the transfromation is legit: res should be zero
%res_1 = X' - (H \ W'*Zica + repmat(mu,1,size(Zica,2)));
%res_2 = Zica -  inv(H\W')*(X' - repmat( mu,1, size(X',2) ) );

%get kymograph:
%Vkymo = transpose(inv(H\W')*V(:,1:embedded_dimension)') ;
%[Kymo_curvdata, Kymo_theta] = Get_Kymograph(Vkymo,8000,eigenvectors,K_delay,K_mode,embedded_dimension);
%Wanring: our modes are randomly sorted every time we run it
len = floor(size(Zica, 1) / 2);

figure
hold on
title('ICA by WenLab')
subplot(1, 2, 1); plot(Zica(2, end - len:end), Zica(5, end - len:end));
xlabel('Z_2')
ylabel('Z_5')
subplot(1, 2, 2); plot(Zica(5, end - len:end), Zica(6, end - len:end));
xlabel('Z_5')
ylabel('Z_6')
hold off
saveas(gcf, fullfile(savefolder,'ICA_phase_2D_by_WenLab.jpg'))

figure
hold on
plot3(-Zica(4, end - len:end), Zica(3, end - len:end), -Zica(1, end - len:end));
xlabel('Z_4')
ylabel('Z_3')
zlabel('Z_1')
hold off
saveas(gcf,fullfile(savefolder,'ICA_phase_3D_by_WenLab.jpg'))
%animation3(X(startf:endf,1),X(startf:endf,2),X(startf:endf,3),1);
%figure; animation3(Zcw(1,startf:endf),Zcw(2,startf:endf),Zcw(3,startf:endf),1);
%figure; animation3(Zica(1,startf:endf),Zica(2,startf:endf),Zica(3,startf:endf),1);
toc


%Estimate the Jacobians and Lyapunov Spectrum
[merr,Jac_data]=learn_jac(XX');
figure
hold on
plot(Jac_data.lexp_avg*ts)
xlabel('Jac-data')
ylabel('lexp*ts')
hold off
saveas(gcf, fullfile(savefolder,'Jacobians.jpg'))


% % Ger periodic orbits and plot all orbits of period 1
% [pertab,pers]=perorb_find_adapt(XX');
% figure
% hold on
% clf;tpper=pertab(pertab(:,6)==1,:);
% plot_clusters_traj(XX,tpper,[3 4 2;6 7 2;1 2 5],ones(12,1),[1 1 1]*0.2)
% xlabel('')
% ylabel('')
% 
% hold off


