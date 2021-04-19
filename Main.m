% User Input Image:
img = imread('StackNinja1.bmp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          Masking Program                          %%%                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nu = img(:,:,2); % Extract the Green Layer.
nu = medfilt2(adapthisteq(nu)); % Increase the contrast adaptively.

% Blur the image to make the image flat.
blur_nu = medfilt2(imgaussfilt(imgaussfilt(nu)));

% Sharpen the image for edge detection.
blur_nu = unsharp_filter(blur_nu);
blur_nu = unsharp_filter(blur_nu);

% padding before edge detection
blur_nu_pad = padarray(blur_nu,[1 1],1,'both');


% Detect the edge.
[Gmag, Gdir] = imgradient(blur_nu_pad,'prewitt');

% Fill the area surrounded by edges, remove the padding as well.
clear_nu = imfill(uint8(Gmag));
clear_nu([1 end],:)=[]; clear_nu(:,[1 end])=[];
Gmag([1 end],:)=[]; Gmag(:,[1 end])=[];
figure;
imshow(clear_nu); 
title('After Detecting the Edge and Filling the Center');

% Plain_nu is plain nuclei, which can be used to fill the hole in edge
% based clear_nu.
plain_nu = imgaussfilt(img(:,:,2));

% Make a outline, this lines is to cut the connection between nuclei.
lines = imfilter(clear_nu,fspecial('log'));
lines3 = imfilter(plain_nu,fspecial('log'))*9 - 40; 


% Merge clear_nu and plain_nu.
figure;
final_nu = clear_nu + plain_nu*0.6;
imshow(final_nu);
title('Final = Clear + Plain');

% Cut Lines from the final nuclei.
final_nu = imfill(final_nu-imgaussfilt(lines3) ,'holes');
figure; imshow(final_nu);
title('Cut the outline between nuclei');

% Open after Cut lines
final_nu = imopen(final_nu,strel('disk',2));
figure; imshow(final_nu);
title('Openning after cut');

% Remove Gray and Binarize
T = showrange(final_nu,60,255);
BW_nu = imbinarize(T,adaptthresh(final_nu));
figure; imshow(BW_nu);
title('Remove grey and binarize adaptively');

BW_nu = imopen(BW_nu,strel('disk',3));
figure; imshow(BW_nu);
title('Masking after openning.');

%bwperim to have the boundary and use blue light to highlight
circles = bwperim(BW_nu);
final_img = img;
final_img(:,:,3) = im2uint8(circles);
figure;imshow(final_img);
title("Boundaries");

figure;imshow(imoverlay(img,BW_nu,'cyan'));
title("Mask Overlay");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       Analysis Program                            %%%                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;imshow(final_img);
CC = bwconncomp(BW_nu);
bbox = regionprops(CC, 'BoundingBox');
rounds = zeros(1,CC.NumObjects);

% Average Brightness Distribution
b_nu = img(:,:,2);
meanBrightnesses = zeros(1,CC.NumObjects);
for ith_nuclei =[1:CC.NumObjects]
    nu_size = numel(CC.PixelIdxList{1,ith_nuclei});
    sumBright = sum(b_nu(CC.PixelIdxList{1,ith_nuclei}));
    meanBright = sumBright/nu_size;
    meanBrightnesses(:,ith_nuclei) = meanBright;
end

for i = 1:CC.NumObjects
    bb=bbox(i).BoundingBox;
    rectangle('Position', bb, 'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 1);
    
    col = bb(1); row = bb(2);
    width = bb(3); height = bb(4);
    
    round = width/height;
    rounds(i) = round;
    
    size = numel(CC.PixelIdxList{i});
    
    text_size = text(col+1,row+1,"Size:" + num2str(size));
    text_bright = text(col+1,row+3,"Brightness:" + num2str(meanBrightnesses(i)));
    text_round = text(col+1, row+5, "Ratio:"+num2str(round));
    
    set(text_size,'Color','w','FontSize', 7, 'FontWeight','bold');
    set(text_bright,'Color','w','FontSize', 7, 'FontWeight','bold');
    set(text_round,'Color','w','FontSize', 7, 'FontWeight','bold');
end
% Size Distribution.
figure;
subplot(1,3,1);
histfit(cellfun(@numel,CC.PixelIdxList),40);
title("Distribution of Nucleus Size");
xlabel('Size of a Nucleus');
ylabel('Nuclei Counts');
sgtitle(sprintf('Nuclei Counts:%d', CC.NumObjects));

subplot(1,3,2);
histfit(meanBrightnesses,40);
title("Distribution of Nucleus Average Brightness");
xlabel('Average Brightness of a Nucleus');
ylabel('Nuclei Counts');

subplot(1,3,3);
histfit(rounds,40);
title("Distribution of Width/Height Ratio");
xlabel('Width/Height Ratio of a Nucleus');
ylabel('Nuclei Counts');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                             Functions                             %%%                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [img] = showrange(x, low, high)
    x(x<low) = 0;
    x(x>high) = 0;
    img = x
end

% Unsharp Filters to sharpen the edge of image.
function [img] = unsharp_filter(x)
    im_blur = imfilter(x, fspecial('gaussian', 10, 20));
    im_edge = x - im_blur;
    img = x + im_edge;
end
