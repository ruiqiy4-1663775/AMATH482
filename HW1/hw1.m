%%
clear; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
average = zeros(64, 64, 64);
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
for j=1:20
   Un(:,:,:)=reshape(Undata(j,:),n,n,n);
   %close all, isosurface(X,Y,Z,abs(Un),0.4)
   %axis([-20 20 -20 20 -20 20]), grid on, drawnow
   ft = fftn(Un);
   average = average + ft;
end
average = abs(fftshift(average)) / max(max(max(abs(average))));
close all, isosurface(Kx,Ky,Kz, average, 0.9)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
[M, I] = max(average(:));


kx0 = Kx(I);
ky0 = Ky(I);
kz0 = Kz(I);

%%
% Define Gaussian Filter

filter = exp(-((Kx - kx0).^2 + (Ky - ky0).^2 + (Kz - kz0).^2));
px = zeros(20);
py = zeros(20);
pz = zeros(20);

for i = 1 : 20
    Un(:, :, :) = reshape(Undata(i, :), n, n, n);
    Un_fft = fftn(Un);
    fftFilter = fftshift(Un_fft) .* filter;
    fftFilter = ifftshift(fftFilter);
    position = ifftn(fftFilter);
    position = position / max(max(max(position)));
    [M, I] = max(position(:));
    px(i) = X(I);
    py(i) = Y(I);
    pz(i) = Z(I);
    isosurface(X,Y,Z,abs(position),0.4)
    axis([-20 20 -20 20 -20 20]), grid on, drawnow
    hold on;
end
xlabel('x coordinates');
ylabel('y coordinates');
zlabel('z coordinates');
title("marble's position");
%%
close all;
plot3(px, py, pz);
xlabel('x coordinates');
ylabel('y coordinates');
zlabel('z coordinates');
title("marble's path");
px(20)
py(20)
pz(20)


