img = imread("image.jpg");
gray = rgb2gray(img);
resized_gray = imresize(gray, [5000, 5000]);
writematrix(resized_gray, 'input.txt', 'Delimiter', 'tab');