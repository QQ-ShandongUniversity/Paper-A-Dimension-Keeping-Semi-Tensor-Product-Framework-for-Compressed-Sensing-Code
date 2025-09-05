% 压缩感知与均值图像MAE分析（归一化版本）
clear; clc; close all;

% 1. 读取原始图像并转换为灰度图
im_original = rgb2gray(imread('image4.bmp')); % 使用您提供的图像
im_original = im2double(im_original); % 将像素值归一化到0-1范围
d_full = size(im_original, 1);

% 设置参数
blockSize = 64; % 图像块大小
numBlocks = 5;  % 随机截取的图像块数量
compressionRatios = 0.05:0.05:0.50; % 压缩率从5%到50%

% 随机选择5个64x64图像块的位置
rng(42); % 设置随机种子以确保结果可重现
blockPositions = zeros(numBlocks, 2);
for i = 1:numBlocks
    row = randi([1, d_full - blockSize + 1]);
    col = randi([1, d_full - blockSize + 1]);
    blockPositions(i, :) = [row, col];
end

% 预分配MAE结果存储
maeResultsReconstructed = zeros(length(compressionRatios), 1);
maeResultsMeanImage = zeros(length(compressionRatios), 1);

% 对每个压缩率进行处理
for crIdx = 1:length(compressionRatios)
    r = compressionRatios(crIdx); % 当前压缩率
    maeReconstructed = zeros(numBlocks, 1);
    maeMeanImage = zeros(numBlocks, 1);
    
    for blockIdx = 1:numBlocks
        % 提取图像块
        row = blockPositions(blockIdx, 1);
        col = blockPositions(blockIdx, 2);
        im_block = im_original(row:row+blockSize-1, col:col+blockSize-1);
        
        % 将图像块向量化（像素值已在0-1范围内）
        x = im_block(:);
        n = length(x);
        
        % 计算均值图像（相邻两个元素求和取平均）
        x_mean = zeros(n/2, 1);
        for i = 1:2:n
            idx = ceil(i/2);
            x_mean(idx) = (x(i) + x(i+1)) / 2;
        end
        
        % 将均值图像扩展回原始尺寸（重复每个值两次）
        x_mean_expanded = repelem(x_mean, 2);
        
        % 设置参数（基于您提供的代码）
        p = 2;
        nf = n/p;
        m = floor(n * r); % 采样数
        
        % 创建测量矩阵（基于您提供的代码）
        row_vector = repmat(p^(-0.5), 1, p);
        phi = randn(m, nf); % 高斯随机矩阵
        phi1 = kron(phi, row_vector);
        
        % 进行压缩感知测量
        y = phi1 * x;
        
        % 创建DCT稀疏基矩阵
        psi = zeros(n);
        for k = 1:n
            t = zeros(n, 1);
            t(k) = 1;
            t = idct(t);
            psi(:, k) = t;
        end
        
        % 组合矩阵
        theta = phi1 * psi;
        
        % 使用l1eq_pd进行重建
        s2 = theta' * y;
        s1 = l1eq_pd(s2, theta, [], y, 1e-3, 20);
        x_reconstructed = psi * s1;
        
        % 确保重建值在0-1范围内
        x_reconstructed = max(0, min(1, x_reconstructed));
        
        % 计算重构图像的MAE
        maeReconstructed(blockIdx) = mean(abs(x_mean_expanded - x_reconstructed));
        
    end
    
    % 存储平均MAE
    maeResultsReconstructed(crIdx) = mean(maeReconstructed);
    
    % 显示进度
    fprintf('完成压缩率 %.0f%% 的处理\n', r*100);
end

% 3. 绘制结果
figure('Position', [100, 100, 1200, 800]);

% 绘制压缩率与MAE关系
subplot(2, 2, 1);
plot(compressionRatios*100, maeResultsReconstructed, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on;
xlabel('Compression ratio (%)');
ylabel('Average MAE');
title('');
grid on;
ylim([0, 0.5]); % 设置Y轴范围，使MAE值在0-1范围内更清晰可见

% 显示原始图像和随机选择的块
subplot(2, 2, 3);
imshow(im_original);
title('');
hold on;
for i = 1:numBlocks
    row = blockPositions(i, 1);
    col = blockPositions(i, 2);
    rectangle('Position', [col, row, blockSize, blockSize], ...
              'EdgeColor', 'r', 'LineWidth', 2);
end
hold off;

