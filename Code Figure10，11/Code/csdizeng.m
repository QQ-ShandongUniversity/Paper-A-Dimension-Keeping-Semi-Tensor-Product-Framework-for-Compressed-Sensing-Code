PSNRCS = zeros(1,10);
%times = zeros(1, 10);
for i = 1:10
im = rgb2gray(imread('image2.bmp'));
%im = imnoise(im,'salt & pepper')
%im =  imnoise(im,'gaussian',0,0.001);
d = size(im,1);
% im = imread('cameraman.tif');
x = double(im(:));
n = length(x);
r = 0.05+0.05*(i-1);
m = floor(n*r);

phi = randn(m,n);
%phi = ToeplitzMtx( m,n ) %toeplitz

%phi = SparseRandomMtx( m,n,16 );

%phi= randi([0,1],m,n);
%phi(phi==0) = -1;
%phi = phi/sqrt(m);%BENOULI
y = phi*x;

psi = zeros(n);
for k = 1:n
    t = zeros(n,1);
    t(k) = 1;
    t = idct(t);
    psi(:,k) = t;
end

theta = phi*psi;


% theta = zeros(m,n);
% 
% for k = 1:n
%     t = zeros(n,1);
%     t(k) = 1;
%     t = idct(t);
%     p = phi*t;
%     theta(:,k) = p;
% end

s2 = theta'*y;

s1 = l1eq_pd(s2, theta, [], y, 1e-3,25);
xp = psi*s1;
%times(i) = toc;
%distance = norm(x - xp);
imRec = reshape(xp,d,d);

%subplot(1,2,1)
%imshow(im,[]);
%xlabel('original');
%subplot(1,2,2)
%imshow(imRec,[]);
%xlabel('reconstructed')

imRec1=uint8(imRec)
PSNR = psnr(im ,imRec1);
PSNRCS(i) = PSNR;
save('CSdizengimage2.mat', 'PSNRCS')
end

