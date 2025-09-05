PSNRDKSTPCS = zeros(1,10);
%times = zeros(1,10);%记录每次运行时间
for i = 1:10;
 im = rgb2gray(imread('image2.bmp'));
%im = imnoise(im,'salt & pepper')
%im = imread('cameraman.tif');
%im =  imnoise(im,'gaussian',0,0.001);
d = size(im,1);
x = double(im(:));%将原始图化为向量
n = length(x);%原始向量尺寸
p = 2;%
nf = n/p;
r = 0.05+0.05*(i-1);%采样率
m = floor(n*r);%采样数
row_vector = repmat(p^(-0.5), 1, p);
%phi = randn(m,nf);%观察矩阵，为高斯随机分布

%phi = ToeplitzMtx( m,nf ) %toeplitz

%phi = SparseRandomMtx( m,nf,16 );

%phi= randi([0,1],m,nf);
%phi(phi==0) = -1;
%phi = phi/sqrt(m);%BENOULi
phi = chaotic_matrix(m, nf);%
phi1 = kron(phi,row_vector);
y = phi1*x;%观察结果

% 稀疏基矩阵，使x=psi*s
psi = zeros(n);%分配空间
for k = 1:n
    t = zeros(n,1);
    t(k) = 1;
    t = idct(t);
    psi(:,k) = t;
end

theta = phi1*psi;%y= theta*s,s是k稀疏的

 %单行避免大矩阵 
 %theta = zeros(m,n);
 
%for k = 1:n
%    t = zeros(n,1);
%     t(k) = 1;
%     t = idct(t);
%     p = phi1*t;
%     theta(:,k) = p;
%end

%初始猜测
s2 = theta'*y;
%tic;%开始计时
s1 = l1eq_pd(s2, theta, [], y, 1e-3,25);%解l1范数
xp = psi*s1;
%time(i) = toc;%存储运行时间
%distance = norm(x - xp);
imRec = reshape(xp,d,d);

%subplot(1,2,1)
%imshow(im,[]);
%xlabel('original');
%subplot(1,2,2);
%imshow(imRec,[]);
%xlabel('reconstructed');
imRec1=uint8(imRec);
PSNR = psnr(im ,imRec1);
PSNRDKSTPCS(i) = PSNR;Toeplitz
save('DKSTPCSdizengChaoticimage2.mat', 'PSNRDKSTPCS')
end
