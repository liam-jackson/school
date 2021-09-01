%Liam Jackson, HW8

clear all, clc;

%3d 
N = 100;
f = zeros(N);
for i = 1:N 
    f(i,:) = .5*(square(2*pi*.05*(1:N), 50)+1);
end

vert = f;
horz = f';

figure(1)
subplot(2,2,1)
%vertical stripes
imshow(vert);

subplot(2,2,2)
%FFT vert
imshow(abs(fft2(vert)));

subplot(2,2,3)
%horizontal stripes
imshow(horz);

subplot(2,2,4)
%FFT horz
imshow(abs(fft2(horz)));

screen = vert .* horz;

figure(2)
imshow(abs(fft2(screen)))


%4a
mouse = imread("labmouse.jpg");
M = double(mouse);

figure(3)
subplot(2,2,1)
imshow("labmouse.jpg");

subplot(2,2,2)
FTM = fft2(M);
FTMcentered = fftshift(FTM);
noise = imshow(FTMcentered,[2500 3500]);

subplot(2,2,3)
mask = zeros(size(FTM)); %I was trying to remove the regions outside the center with a rotated "T" shape. It didn't really work at all but I ran out of time to troubleshoot. 
for i = 240:242
    for j = 1:100
        mask(i,j) = 1;
    end
end
for k = 1:482
    for l = 143:144
        mask(k,l) = 1;
    end
end
filtered = FTM .* mask; 
post = ifft2(filtered);
imshow(post,[])

subplot(2,2,4)
FTpostcent = fftshift(fft2(post));
imshow(FTpostcent)

%5

figure(4)
f = imread('toy-story.jpg');
f = rgb2gray(f);
f = double(f);

F = fftshift(fft2(f));

subplot(2,2,1)
imshow(f,[]);

subplot(2,2,2)
imshow(abs(F))
caxis([0 numel(f)])

subplot(2,2,3)
sorted_values = sort(abs(F(:)));
threshold = sorted_values(round(.95*end));
Fcompressed = F;
Fcompressed(abs(F)<threshold) = 0;

fcompressed = abs(ifft2(Fcompressed));
imshow(fcompressed,[])

subplot(2,2,4)
imshow(fftshift(fft2(fcompressed)),[]);


figure(5)
f = imread('monsters-inc.jpg');
f = rgb2gray(f);
f = double(f);

F = fftshift(fft2(f));

subplot(2,2,1)
imshow(f,[]);

subplot(2,2,2)
imshow(abs(F))
caxis([0 numel(f)])

subplot(2,2,3)
sorted_values = sort(abs(F(:)));
threshold = sorted_values(round(.95*end));
Fcompressed = F;
Fcompressed(abs(F)<threshold) = 0;

fcompressed = abs(ifft2(Fcompressed));
imshow(fcompressed,[])

subplot(2,2,4)
imshow(fftshift(fft2(fcompressed)),[]);











