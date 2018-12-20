mex bestexemplarhelper.c
% [i1,i2,i3,c,d,mov] = inpaint('dataset/nyu_img_601.png','dataset/nyu_depth_601.png','dataset/nyu_mask_601.png',[0 255 0], 0.5);
% plotall

fileFolder = fullfile('dataset');
imgFiles = dir(fullfile(fileFolder,'nyu_img_*.png'));
filename = {imgFiles.name}';
for i = 1:length(filename)
    img = strcat(fileFolder,'/', filename{i});
    depth = strrep(img, 'img', 'depth');
    mask = strrep(img, 'img', 'mask');
    inpaint(img, depth, mask, [0 255 0], 0.4);
end