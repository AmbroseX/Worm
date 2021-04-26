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