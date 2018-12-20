function [inpaintedImg,origImg,fillImg,C,D,fillMovie] = inpaint(imgFilename,depthFilename,fillFilename,fillColor,scaleFactor)
%INPAINT  Exemplar-based inpainting.
%
% Usage:   [inpaintedImg,origImg,fillImg,C,D,fillMovie] ...
%                = inpaint(imgFilename,fillFilename,fillColor)
% Inputs: 
%   imgFilename    Filename of the original image.
%   fillFilename   Filename of the image specifying the fill region. 
%   fillColor      1x3 RGB vector specifying the color used to specify
%                  the fill region.
% Outputs:
%   inpaintedImg   The inpainted image; an MxNx3 matrix of doubles. 
%   origImg        The original image; an MxNx3 matrix of doubles.
%   fillImg        The fill region image; an MxNx3 matrix of doubles.
%   C              MxN matrix of confidence values accumulated over all iterations.
%   D              MxN matrix of data term values accumulated over all iterations.
%   fillMovie      A Matlab movie struct depicting the fill region over time. 
%
% Example:
%   [i1,i2,i3,c,d,mov] = inpaint('img.png','depth.png','mask.png',[0 255 0],scale);
%   plotall;           % quick and dirty plotting script
%   close; movie(mov); % grab some popcorn 
%
%   author: Sooraj Bhat
%   Modified by Marcel Davey & John Gu on
%   11/30/05 to run on Matlab 7.0.4.365 (R14) Service Pack 2
   

warning off MATLAB:divideByZero
[img_hd,depth_hd,fillImg_hd,fillRegion_hd] = loadimgs(imgFilename,depthFilename,fillFilename,fillColor,scaleFactor);

img = imresize(img_hd, scaleFactor);
depth = imresize(depth_hd, scaleFactor);
fillImg = imresize(fillImg_hd, scaleFactor);
fillRegion = imresize(fillRegion_hd, scaleFactor);


img = double(img);
depth = double(depth);
origImg = img;
origDepth = depth;
origFillRegion = fillRegion;
H_hd = size(origImg, 1);
W_hd = size(origImg, 2);

% Gradient map
[Dx Dy] = gradient(depth);
depth_gradient = cat(3, Dx, Dy);
origGradient = depth_gradient;

ind = img2ind(img);
ind_hd = img2ind(img_hd);        % High resolution index map
sz = [size(img,1) size(img,2)];
sourceRegion = ~fillRegion;


% Initialize isophote values
[Ix(:,:,3) Iy(:,:,3)] = gradient(img(:,:,3));
[Ix(:,:,2) Iy(:,:,2)] = gradient(img(:,:,2));
[Ix(:,:,1) Iy(:,:,1)] = gradient(img(:,:,1));
Ix = sum(Ix,3)/(3*255); Iy = sum(Iy,3)/(3*255);
temp = Ix; Ix = -Iy; Iy = temp;  % Rotate gradient 90 degrees to get isophote

IDx = -Dy; IDy = Dx; % Isophote of depth

% Initialize confidence and data terms
C = double(sourceRegion);
D = repmat(-.1,sz);
iter = 1;
% Visualization stuff
if nargout==6
  fillMovie(1).cdata=uint8(img); 
  fillMovie(1).colormap=[];
  origImg(1,1,:) = fillColor;
  iter = 2;
end

% Seed 'rand' for reproducible results (good for testing)
rand('state',0);
% Parameter for priority: lambda * RGB + (1 - lambda) * depth
lambda = 0.8;

% Loop until entire fill region has been covered
while any(fillRegion(:))
  % Find contour & normalized gradients of fill region
  fillRegionD = double(fillRegion); % Marcel 11/30/05
  dR = find(conv2(fillRegionD,[1,1,1;1,-8,1;1,1,1],'same')>0); % Marcel 11/30/05
 %dR = find(conv2(fillRegion,[1,1,1;1,-8,1;1,1,1],'same')>0);  % Original
  
  [Nx,Ny] = gradient(double(~fillRegion)); % Marcel 11/30/05
 %[Nx,Ny] = gradient(~fillRegion);         % Original
  N = [Nx(dR(:)) Ny(dR(:))];
  N = normr(N);  
  N(~isfinite(N))=0; % handle NaN and Inf
  
  % Compute confidences along the fill front
  for k=dR'
    Hp = getpatch(sz,k);
    q = Hp(~(fillRegion(Hp)));  % Get source region in the patch
    C(k) = sum(C(q))/numel(Hp);
  end
  
  % Compute patch priorities = confidence term * data term
  % D(dR) = abs(Ix(dR).*N(:,1)+Iy(dR).*N(:,2)) + 0.001;      % Added a little constant to data term.
  D(dR) = lambda * abs(Ix(dR).*N(:,1)+Iy(dR).*N(:,2)) + (1-lambda) * abs(IDx(dR).*N(:,1)+IDy(dR).*N(:,2)) + 0.001;        % Modified data term for depth. 
  priorities = C(dR).* D(dR);
  
  % Find patch with maximum priority, Hp
  [unused,ndx] = max(priorities(:));
  p = dR(ndx(1));
  [Hp,rows,cols] = getpatch(sz,p);
  toFill = fillRegion(Hp);
  
  % Find exemplar that minimizes error, Hq
  [Hq, rowsq, colsq] = bestexemplar(img,img(rows,cols,:),depth_gradient, depth_gradient(rows,cols,:),toFill',sourceRegion, rows(1), cols(1));
%  [Hq, rowsq, colsq] = bestexemplar(img,img(rows,cols,:),depth, depth(rows,cols),toFill',sourceRegion);

  % [Hq, rowsq, colsq] = bestexemplar_depth(img, depth, img(rows,cols,:), depth(rows,cols), toFill', sourceRegion);
  
  % Update fill region
  toFill = logical(toFill);                 % Marcel 11/30/05
  fillRegion(Hp(toFill)) = false;
  
  % Propagate confidence & isophote values
  C(Hp(toFill))  = C(p);
  Ix(Hp(toFill)) = Ix(Hq(toFill));
  Iy(Hp(toFill)) = Iy(Hq(toFill));
  
  % Copy image data from Hq to Hp
  ind(Hp(toFill)) = ind(Hq(toFill));
  img(rows,cols,:) = ind2img(ind(rows,cols),origImg);
  depth(rows,cols) = origDepth(ind(rows,cols));
  depth_gradient(rows,cols,:) = ind2img_n(ind(rows,cols),origGradient,2);
  
  % high resolution
  rows_hd = floor(rows(1)/scaleFactor)-1:floor(rows(1)/scaleFactor)+floor((length(rows)-1)/scaleFactor);
  cols_hd = floor(cols(1)/scaleFactor)-1:floor(cols(1)/scaleFactor)+floor((length(cols)-1)/scaleFactor);
  rowsq_hd = floor(rowsq(1)/scaleFactor)-1:floor(rowsq(1)/scaleFactor)+floor((length(rowsq)-1)/scaleFactor);
  colsq_hd = floor(colsq(1)/scaleFactor)-1:floor(colsq(1)/scaleFactor)+floor((length(colsq)-1)/scaleFactor);
  ind_hd(rows_hd, cols_hd) = ind_hd(rowsq_hd, colsq_hd);
  
  % Visualization stuff
  if nargout==6
    ind2 = ind;
    ind2(logical(fillRegion)) = 1;          % Marcel 11/30/05
    %ind2(fillRegion) = 1;                  % Original
    fillMovie(iter).cdata=uint8(ind2img(ind2,origImg)); 
    fillMovie(iter).colormap=[];
  end
  iter = iter+1;
end

inpaintedImg=img;

% Depth
reconstructedDepth = reconstruct(depth, origFillRegion, depth_gradient(:,:,1), depth_gradient(:,:,2));
reconstructedDepthScaled = imresize(reconstructedDepth, 'OutputSize', size(depth_hd));
fillDepth = depth_hd;
fillDepth(fillRegion_hd) = reconstructedDepthScaled(fillRegion_hd);

imwrite(uint8(fillDepth), strcat('inpainted_', depthFilename));
% High resolution
inpaintedImgHD=ind2img(ind_hd, img_hd);

imwrite(uint8(inpaintedImgHD), strcat('inpainted_', imgFilename));
% imwrite(uint8(inpaintedImgHD), 'output_hd.png');

% imshow(depth_gradient(:,:,1));

%---------------------------------------------------------------------
% Scans over the entire image (with a sliding window)
% for the exemplar with the lowest error. Calls a MEX function.
%---------------------------------------------------------------------
function [Hq rows cols] = bestexemplar(img,Ip,depth,Dp,toFill,sourceRegion,center_y,center_x)
m=size(Ip,1); mm=size(img,1); n=size(Ip,2); nn=size(img,2);
best = bestexemplarhelper(mm,nn,m,n,img,Ip,depth,Dp,toFill,sourceRegion,center_y, center_x);
Hq = sub2ndx(best(1):best(2),(best(3):best(4))',mm);
rows = best(1):best(2);
cols = best(3):best(4);

% New bestexemplar

function [Hq, rows, cols] = bestexemplar_depth(img, depth, Ip, Dp, toFill, sourceRegion)
imageSize = size(img);
depthSize = size(depth);
assert(~any(depthSize(1:2) ~= imageSize(1:2)),'Depth map should have the same size with image');

H = imageSize(1);
W = imageSize(2);
h = size(Ip, 1);         % size of patch
w = size(Ip, 2);

Filled = ~toFill;
alpha = 1.54e-4;           % Scalar parameter for depth

errorMap = Inf(H-h+1, W-w+1);

% for all patches
parfor j = 1:H-h+1
    for i = 1:W-w+1
        if (any(any(~sourceRegion(j:j+h-1,i:i+w-1))))  % Skip Fillregion
            continue; 
        end
        templatePatch = img(j:j+h-1,i:i+w-1,:);
        templatePatchDepth = depth(j:j+h-1,i:i+w-1);
        errorMap(j,i) = sum(sum(((templatePatch(:,:,1) - Ip(:,:,1)).*Filled).^2)) ...
                      + sum(sum(((templatePatch(:,:,2) - Ip(:,:,2)).*Filled).^2)) ...
                      + sum(sum(((templatePatch(:,:,3) - Ip(:,:,3)).*Filled).^2)) ...
              + alpha * sum(sum(((templatePatchDepth - Dp).*Filled).^2));        % calculate error (L2 norm over all filled pix within the patch)
    end
end

[argvalue, argmin] = min(errorMap(:));
[J,I] = ind2sub(size(errorMap), argmin);

rows = J:J+h-1;
cols = I:I+w-1;
Hq = sub2ndx(J:J+h-1,(I:I+w-1)',H);

%---------------------------------------------------------------------
% Returns the indices for a 9x9 patch centered at pixel p.
%---------------------------------------------------------------------
function [Hp,rows,cols] = getpatch(sz,p)
% [x,y] = ind2sub(sz,p);  % 2*w+1 == the patch size
w=4; p=p-1; y=floor(p/sz(1))+1; p=rem(p,sz(1)); x=floor(p)+1;
rows = max(x-w,1):min(x+w,sz(1));
cols = (max(y-w,1):min(y+w,sz(2)))';
Hp = sub2ndx(rows,cols,sz(1));


%---------------------------------------------------------------------
% Converts the (rows,cols) subscript-style indices to Matlab index-style
% indices.  Unforunately, 'sub2ind' cannot be used for this.
%---------------------------------------------------------------------
function N = sub2ndx(rows,cols,nTotalRows)
X = rows(ones(length(cols),1),:);
Y = cols(:,ones(1,length(rows)));
N = X+(Y-1)*nTotalRows;


%---------------------------------------------------------------------
% Converts an indexed image into an RGB image, using 'img' as a colormap
%---------------------------------------------------------------------
function img2 = ind2img(ind,img)
for i=3:-1:1, temp=img(:,:,i); img2(:,:,i)=temp(ind); end;

function img2 = ind2img_n(ind,img,channel)
for i=channel:-1:1, temp=img(:,:,i); img2(:,:,i)=temp(ind); end;
%---------------------------------------------------------------------
% Converts an RGB image into a indexed image, using the image itself as
% the colormap.
%---------------------------------------------------------------------
function ind = img2ind(img)
s=size(img); ind=reshape(1:s(1)*s(2),s(1),s(2));


%---------------------------------------------------------------------
% Loads the an image and it's fill region, using 'fillColor' as a marker
% value for knowing which pixels are to be filled.
%---------------------------------------------------------------------
function [img,depth,fillImg,fillRegion] = loadimgs(imgFilename,depthFilename,fillFilename,fillColor,scaleFactor)
img = imread(imgFilename); fillImg = imread(fillFilename);
depth = imread(depthFilename);

% NYU Depth crop

img = img(7:474, 8:632, :);
fillImg = fillImg(7:474, 8:632, :);
depth = depth(7:474, 8:632);

if(size(fillImg,3)>1)
    fillRegion = fillImg(:,:,1)==fillColor(1) & fillImg(:,:,2)==fillColor(2) & fillImg(:,:,3)==fillColor(3);
else
    fillRegion = fillImg;
end

% fillRegion = fillImg(:,:,1)==fillColor(1) & ...
%     fillImg(:,:,2)==fillColor(2) & fillImg(:,:,3)==fillColor(3);

function [A] = normr(N)
    for ii=1:size(N,1)
        A(ii,:) = N(ii,:)/norm(N(ii,:));
    end

