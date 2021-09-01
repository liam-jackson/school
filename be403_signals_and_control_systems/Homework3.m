%Liam Jackson, HW3

%Question 1e

x = [1 2 0 -1]';
h = [1 -1 1 -1]';

yall = conv(h,x);
y = yall([1:4],1)

clear

%Question 3 check

x = [1 0 1; 0 1 0; 1 0 1];
h = [0 .1 0; .1 .6 .1; 0 .1 0];
yall = conv2(x,h);
y = yall([2:4],[2:4])

clear

%Question 4

pic = single(imread("cells.jpg"));
gray = mat2gray(pic);

h1 = [0 -1 0; 
    -1 5 -1; 
    0 -1 0];
h2 = [0 1 0;
    1 -4 1;
    0 1 0];
h3 = [-2 -1 0;
    -1 1 1;
    0 1 2];

%a
figure
subplot(2,2,1)
imshow(gray)
title('Grayscale')

%b
subplot(2,2,2)
pic2 = conv2(h1,gray);
imshow(pic2)
title('Sharpened')

%c
subplot(2,2,3)
pic3 = conv2(h2,gray);
imshow(pic3,[])
title('Edge Detection')

%d
subplot(2,2,4)
pic4 = conv2(h3,gray);
imshow(pic4)
title('Embossed')

%e
%{
Based on some cursory research, this property denotes a Normalized kernel. 
This maintains the average intensity of a pixel across the input signal 
and the generated signal.  
%}



