function out = preproc_PIV (in,roirect,clahe,num_tiles,highp,highpsize,intenscap,wienerwurst,wienerwurstsize,minintens,maxintens)
% preprocess the images

%% initialize
% grayscale, only one channel
if size(in,3)>1
    in(:,:,2:end)=[];
end

% decode roirect into x, y, width, height
% numel: number of elements in an array or subscripted array expression
% x width: horizontal; y height: vertical;
if numel(roirect)>0
    x=roirect(1);
    y=roirect(2);
    width=roirect(3);
    height=roirect(4);
else % default
    x=1;
    y=1;
    width=size(in,2)-1;
    height=size(in,1)-1;
end

% roi: region of interest of the image, (x,y,width,height)
in_roi=in(y:y+height,x:x+width);

% adjust histogram, for 8bit and 16 bit
B = imadjust(in_roi, [minintens;maxintens],[]); % if defaults = 0 and 1, then no effect
B=(double(B)/max(max(double(B))));

in_roi=uint8(B*255); % convert back to uint8

%% intensity capping
if intenscap == 1
    % Intensity Capping: a simple method to improve cross-correlation PIV results
    n = 2; 
    up_lim_im_1 = median(double(in_roi(:))) + n*std2(in_roi); % upper limit for the image
    brightspots_im_1 = in_roi > up_lim_im_1;     % bright spots in the image
    capped_im_1 = in_roi; 
    capped_im_1(brightspots_im_1) = up_lim_im_1; % capped image
    in_roi=capped_im_1;
end

%% CLAHE
if clahe == 1
    in_roi=adapthisteq(in_roi,'NumTiles',num_tiles,...
        'ClipLimit',0.01,'NBins',256,'Range','full','Distribution','uniform');
end

%% highpass
if highp == 1
    h = fspecial('gaussian',highpsize,highpsize);
    in_roi=double(in_roi-(imfilter(in_roi,h,'replicate')));
    in_roi=in_roi/max(max(in_roi))*255;
end

if wienerwurst == 1
    in_roi=wiener2(in_roi,[wienerwurstsize wienerwurstsize]);
end

out=in;
out(y:y+height,x:x+width)=in_roi;
out=uint8(out);
