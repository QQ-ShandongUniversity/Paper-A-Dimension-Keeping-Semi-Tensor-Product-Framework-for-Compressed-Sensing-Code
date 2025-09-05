% 图像处理与误差分析脚本
clear; clc; close all;

% 1. 读取原始图像并转换为灰度图（如果不是灰度图）
originalImage = imread('image4.bmp'); % 请替换为您的图像路径
if size(originalImage, 3) == 3
    originalImage = rgb2gray(originalImage);
end

% 将图像转换为双精度类型以便计算
originalImage = im2double(originalImage);

% 2. 将图像转换为向量信号
imageVector = originalImage(:); % 将图像展平为列向量

% 3. 对相邻两个灰度值求平均并重构信号
% 确保向量长度为偶数
if mod(length(imageVector), 2) ~= 0
    imageVector = imageVector(1:end-1); % 如果是奇数，去掉最后一个元素
end

% 将向量重塑为2行矩阵，便于计算每对元素的平均值
pairedVector = reshape(imageVector, 2, []);

% 计算每对相邻像素的平均值
averagedValues = mean(pairedVector, 1);

% 重构信号：将每个平均值重复两次
reconstructedVector = repelem(averagedValues, 2)';

% 4. 计算MAE值
maeValue = mean(abs(imageVector - reconstructedVector));
fprintf('平均绝对误差(MAE): %.4f\n', maeValue);

% 5. 计算每个像素的绝对误差
pixelErrors = abs(imageVector - reconstructedVector);

% 6. 绘制误差分布直方图
figure;
subplot(2, 2, 1);
imshow(originalImage);
title('Original image');

subplot(2, 2, 2);
% 将重构向量重新转换为图像格式以便显示
reconstructedImage = reshape(reconstructedVector, size(originalImage, 1), size(originalImage, 2));
imshow(reconstructedImage);
title('Reconstructed image');

subplot(2, 2, 3);
errorImage = reshape(pixelErrors, size(originalImage, 1), size(originalImage, 2));
imshow(errorImage, []);
colormap('jet');
colorbar;
title('Heat map of error distribution');

subplot(2, 2, 4);
histogram(pixelErrors, 50); % 使用50个bin
xlabel('Absolute error value');
ylabel('Frequency');
title('Error distribution histogram');
grid on;

% 添加统计信息到直方图
meanError = mean(pixelErrors);
medianError = median(pixelErrors);
maxError = max(pixelErrors);

text(0.7*max(pixelErrors), 0.8*max(ylim), sprintf('MAE: %.4f\nmedian value: %.4f\nmaximum value: %.4f', ...
    meanError, medianError, maxError), 'VerticalAlignment', 'top');

% 调整布局
set(gcf, 'Position', [100, 100, 1200, 800]);