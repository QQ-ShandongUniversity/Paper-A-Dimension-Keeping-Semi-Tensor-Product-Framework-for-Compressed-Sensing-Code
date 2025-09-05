%PSNRSTPCS = zeros(1,50);
%times = zeros(1, 50);%记录每次运行时间
%for i = 1:50;
im = rgb2gray(imread('image4.bmp'));
%im = imnoise(im,'salt & pepper');
% im = imread('cameraman.tif');
im =  imnoise(im,'gaussian',0,0.001);
d = size(im,1);
x = double(im(:));%将原始图化为向量
n = length(x);%原始向量尺寸
p = 2;
nf = n/p;
r = 0.5;%采样率
m = floor(n*r);%采样数
mf = m/p;
phi = randn(mf,nf);%观察矩阵，为高斯随机分布
phi2 = kron(phi,eye(p));
y = phi2*x;%观察结果

% 稀疏基矩阵，使x=psi*s
psi = zeros(n);%分配空间
for k = 1:n
    t = zeros(n,1);
    t(k) = 1;
    t = idct(t);
    psi(:,k) = t;
end

theta = phi2*psi;%y= theta*s,s是k稀疏的

% 单行避免大矩阵 
% theta = zeros(m,n);
% 
% for k = 1:n
%     t = zeros(n,1);
%     t(k) = 1;
%     t = idct(t);
%     p = phi*t;
%     theta(:,k) = p;
% end

%初始猜测
s2 = theta'*y;
%tic;%开始计时
s1 = l1eq_pd(s2, theta, [], y, 1e-3,25);%解l1范数
xp = psi*s1;
%times(i) = toc;%存储运行时间
imRec = reshape(xp,d,d);

subplot(1,2,1)
imshow(im,[]);
xlabel('original');
subplot(1,2,2)
imshow(imRec,[]);
xlabel('reconstructed')

imRec1 = uint8(imRec);
PSNR = psnr(im ,imRec1);
%PSNRSTPCS(i) = PSNR;
%save('STPCSgaussimage4.mat', 'PSNRSTPCS')
%end

