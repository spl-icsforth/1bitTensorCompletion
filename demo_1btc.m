%%%%%%%%%%%%%%%%%%%%    1-Bit Tensor Completion    %%%%%%%%%%%%%%%%%%%%
clear
close all
format long

addpath('.\inexact_alm_mc');
addpath('.\inexact_alm_mc\PROPACK');

rand('seed',2013);

%% real data
% Hyperspectral image - Indian Pines (145x145x200)
load('.\Data\Indian_pines_corrected.mat')
X = indian_pines_corrected;
M = imresize(X,1/2,'nearest'); % 73x73x200

Nway = size(M);
N = ndims(M);

% R for PSNR
bits_orig = ceil(log2(max(M(:))));
R = 2^bits_orig-1;

%% Quantization
bit = 1;

% Different quantization for each pixel
Mq = zeros(size(M));
part = zeros(Nway(1),Nway(2),2^bit-1);
codebook = zeros(Nway(1),Nway(2),2^bit);
mn = zeros(Nway(1),Nway(2));
mx = zeros(Nway(1),Nway(2));
for j = 1:Nway(1)
    for k = 1:Nway(2)
        [Mq(j,k,:),part(j,k,:),codebook(j,k,:),mn(j,k),mx(j,k)] = quantization(M(j,k,:),bit);
    end
end

% Estimate the tensors U and L, the upper and lower bin boundaries
U = zeros(size(M));
L = zeros(size(M));
for j = 1:Nway(1)
    for k = 1:Nway(2)
        [U(j,k,:),L(j,k,:)] = bin_boundaries(Mq(j,k,:),part(j,k,:),codebook(j,k,:),mn(j,k),mx(j,k));
    end
end

%% 
% Percentage of measurments
A = 5:5:100;

% Model: 1) Logistic, 2) Probit
model = 2;

% Tensor of indices = midx
h = 0;
midx = zeros(Nway);
for k = 1:Nway(3)
    for j = 1:Nway(2)
        midx(:,j,k) = (h+1):(h+Nway(1));
        h = h+Nway(1);
    end
end

%% 
Mrec = cell(length(A),1);
Mhat_ten = cell(length(A),N);

nmse = zeros(length(A),N);
nmse_av = zeros(length(A),1);
psnr_er = zeros(length(A),N);
psnr_averr = zeros(length(A),1);

for i = 1:length(A)
    % Observe 'A(i)' percent of the entries in Y
    fprintf('Sampling percentage: %d%%\n', A(i));
    p = round((A(i)/100)*prod(Nway));
    idx = randsample(prod(Nway),p);
    data1 = Mq(idx);
    
    % Measurement tensor
    Y = zeros(size(M(:)));
    Y(idx) = data1;
    Y = reshape(Y,size(M));
    
    % Dynamic weights
    alpha = zeros(N,1);
    fit = zeros(N,1);
    
    for n = 1:N
        fprintf('Sampling percentage: %d%%, Unfolding %d\n', A(i),n);

        % Unfold the tensor of measurements
        Yn = Unfold(Y, Nway, n);
        
        % Unfold the tensor of indices
        midxn = Unfold(midx,Nway,n);
        % Find the indices of observations for each unfolding
        idxn = zeros(length(idx),1);
        for j = 1:length(idx)
            idxn(j) = find(midxn(:) == idx(j));
        end

        % Unfold the tensors of bin boundaries
        Un = Unfold(U, Nway, n);
        Ln = Unfold(L, Nway, n);

        % Recover the estimated matrix Mhat using the quantized matrix completion algorithm QMC
        Mhat = QMC(Yn,idxn,Un,Ln,model);

        % Fold the tensor
        Mhat_ten{i,n} = Fold(Mhat,Nway,n);
        
        % Compute the fitting error and the weight alpha
        fit(n) = norm(Mhat_ten{i,n}(idx)-data1);
        alpha(n) = fit(n)^(-1);

        % Compute error for each unfolding
        Mn = Unfold(M, Nway, n);
        % Normalized mean square error
        nmse(i,n) = norm(Mhat-Mn,'fro')/norm(Mn,'fro');
        % PSNR
        peaksnr = zeros(Nway(3),1);
        for j = 1:Nway(3)
            peaksnr(j) = psnr(Mhat_ten{i,n},M,R);
        end
        psnr_er(i,n) = mean(peaksnr);
    end

    alpha = alpha/sum(alpha);
    % Estimate the weighted sum of the recovered tensors Mhat_ten
    Mrec{i} = zeros(Nway);
    for n = 1:N
        Mrec{i} = Mrec{i}+alpha(n)*Mhat_ten{i,n};
    end
    Mrec{i} = round(Mrec{i});
    
    % Compute error
    % Normalized mean square error
    nmse_av(i) = norm(Mrec{i}(:)-M(:))/norm(M(:));
    % PSNR
    peaksnr = zeros(Nway(3),1);
    for j = 1:Nway(3)
        peaksnr(j) = psnr(Mrec{i},M,R);
    end
    psnr_averr(i) = mean(peaksnr);
end
%% Plot the results: error - number of measurements %%
% Error for each unfolding
figure;
plot(A, nmse(:,1), 'r-*');
hold on
plot(A, nmse(:,2), 'b-*');
hold on
plot(A, nmse(:,3), 'g-*');
xlabel('Sampling percentage')
ylabel('Normalized Mean Square Error')
ylim([0 1])

figure;
plot(A, psnr_er(:,1), 'r-*');
hold on
plot(A, psnr_er(:,2), 'b-*');
hold on
plot(A, psnr_er(:,3), 'g-*');
xlabel('Sampling percentage')
ylabel('PSNR in dB')

% Error of the estimated tensor
figure;
plot(A, nmse_av, '-*');
xlabel('Sampling percentage')
ylabel('Normalized Mean Square Error')
ylim([0 1])

figure;
plot(A, psnr_averr, '-*');
xlabel('Sampling percentage')
ylabel('PSNR in dB')
