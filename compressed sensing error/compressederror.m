% 压缩感知与MAE分析
clear; clc; close all;

% 1. 读取原始图像并转换为灰度图
originalImage = imread('image4.bmp'); % 请替换为您的图像路径
if size(originalImage, 3) == 3
    originalImage = rgb2gray(originalImage);
end
originalImage = im2double(originalImage); % 转换为双精度

% 设置参数
blockSize = 64; % 图像块大小
numBlocks = 5;  % 随机截取的图像块数量
compressionRatios = 0.05:0.05:1.0; % 压缩率从5%到100%
mValues = 0.05:0.05:0.50; % m从5%到50%

% 获取图像尺寸
[height, width] = size(originalImage);

% 随机选择5个64x64图像块的位置
blockPositions = zeros(numBlocks, 2);
for i = 1:numBlocks
    row = randi([1, height - blockSize + 1]);
    col = randi([1, width - blockSize + 1]);
    blockPositions(i, :) = [row, col];
end

% 预分配MAE结果存储
maeResults = zeros(length(compressionRatios), 1);
maeDiffResults = zeros(length(mValues), 1);

% 对每个压缩率进行处理
for crIdx = 1:length(compressionRatios)
    cr = compressionRatios(crIdx);
    m = round(blockSize^2 * cr); % 测量次数
    
    % 对每个图像块进行处理
    blockMae = zeros(numBlocks, 1);
    
    for blockIdx = 1:numBlocks
        % 提取图像块
        row = blockPositions(blockIdx, 1);
        col = blockPositions(blockIdx, 2);
        imageBlock = originalImage(row:row+blockSize-1, col:col+blockSize-1);
        
        % 向量化图像块
        x = imageBlock(:);
        n = length(x);
        
        % 生成高斯随机测量矩阵
        Phi = randn(m, n) / sqrt(m);
        
        % 进行压缩感知测量
        y = Phi * x;
        
        % 使用DCT进行稀疏表示和重建
        Psi = dctmtx(n)'; % DCT变换矩阵
        A = Phi * Psi;    % 组合矩阵
        
        % 使用L1最小化进行重建 (使用MATLAB的l1eq_pd函数，需要安装l1magic工具箱)
        % 如果没有安装l1magic，可以使用其他优化方法，如OMP
        theta_hat = l1eq_pd(zeros(n,1), A, [], y, 1e-3,20);
        
        % 重建信号
        x_hat = Psi * theta_hat;
        
        % 计算MAE
        blockMae(blockIdx) = mean(abs(x - x_hat));
    end
    
    % 存储平均MAE
    maeResults(crIdx) = mean(blockMae);
    
    % 计算MAE差值 (对于m从5%到50%)
    if crIdx <= length(mValues)
        cr2 = 2 * mValues(crIdx);
        if cr2 <= 1
            cr2Idx = find(compressionRatios == cr2);
            if ~isempty(cr2Idx)
                maeDiffResults(crIdx) = maeResults(crIdx) - maeResults(cr2Idx);
            end
        end
    end
end

% 2. 绘制结果
figure('Position', [100, 100, 1200, 800]);

% 绘制压缩率与MAE关系
subplot(2, 2, 1);
plot(compressionRatios*100, maeResults, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
xlabel('Compression ratio (%)');
ylabel('Average MAE');
title('');
grid on;

% 绘制MAE差值
subplot(2, 2, 2);
bar(mValues*100, maeDiffResults, 'FaceColor', [0.5, 0.5, 0.5]);
xlabel('Compression ratio (%)');
ylabel('MAE(m) - MAE(2m)');
title('');
grid on;

% 显示原始图像和随机选择的块
subplot(2, 2, 3);
imshow(originalImage);
title('');
hold on;
for i = 1:numBlocks
    row = blockPositions(i, 1);
    col = blockPositions(i, 2);
    rectangle('Position', [col, row, blockSize, blockSize], ...
              'EdgeColor', 'r', 'LineWidth', 2);
end
hold off;

% 显示一个示例图像块的重建结果 (使用最高压缩率)
subplot(2, 2, 4);
row = blockPositions(1, 1);
col = blockPositions(1, 2);
imageBlock = originalImage(row:row+blockSize-1, col:col+blockSize-1);

% 使用最低压缩率(5%)重建示例块
m_low = round(blockSize^2 * 0.05);
x = imageBlock(:);
n = length(x);
Phi_low = randn(m_low, n) / sqrt(m_low);
y_low = Phi_low * x;
Psi = dctmtx(n)';
A_low = Phi_low * Psi;
theta_hat_low = l1eq_pd(zeros(n,1), A_low, [], y_low, 1e-3);
x_hat_low = Psi * theta_hat_low;
imageBlockLow = reshape(x_hat_low, blockSize, blockSize);

% 使用最高压缩率(100%)重建示例块
m_high = round(blockSize^2 * 1.0);
Phi_high = randn(m_high, n) / sqrt(m_high);
y_high = Phi_high * x;
A_high = Phi_high * Psi;
theta_hat_high = l1eq_pd(zeros(n,1), A_high, [], y_high, 1e-3);
x_hat_high = Psi * theta_hat_high;
imageBlockHigh = reshape(x_hat_high, blockSize, blockSize);

% 显示对比
imshowpair(imageBlockLow, imageBlockHigh, 'montage');
title('左: 5%压缩率重建, 右: 100%压缩率重建');

% 保存结果
save('compression_mae_results.mat', 'compressionRatios', 'maeResults', 'mValues', 'maeDiffResults');