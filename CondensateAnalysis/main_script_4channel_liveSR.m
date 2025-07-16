clc;close all;clear;
%%

filepath = 'D:/ImageData/20250512-pro-jq1-2h/single-cell/';

filename_list = dir([filepath, '*.tif']);

pixelSize = 94; %nm
%%
seg_points = [5, 5];
for file_iter = 1:length(filename_list)

filename = [filename_list(file_iter).name(1:(end-4)), '.tif'];
mkdir([filepath, filename(1:(end-4))]);

z_stack = 1; c_channel = 4;
channel_labels = {'DNA', 'RNA', 'OCT4', 'BRD4'};

% read image series
temp = TIFFreader([filepath, filename]);
h = size(temp, 1); w = size(temp, 2);
numberOfPages = size(temp, 3)/z_stack/c_channel;

img_series = uint16(zeros(h, w, numberOfPages, z_stack, c_channel));
for frame_iter = 1:numberOfPages
    for z_iter = 1:z_stack
        for c_iter = 1:c_channel
            img_series(:, :, frame_iter, z_iter, c_iter) = temp(:, :, (frame_iter-1)*z_stack*c_channel+(z_iter-1)*c_channel+c_iter);
        end
    end
end

% select out max-intensity Z layer
img_series_max = uint16(zeros(h, w, numberOfPages, c_channel));
for frame_iter = 1:numberOfPages
    img_series_max(:, :, frame_iter, :) = img_series(:, :, frame_iter, 1, :);
end
clear img_series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% nucleus mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c_iter = 3; % using OCT4 channel
img_series = img_series_max(:, :, :, c_iter);
% preprocessing parameter
gaussian_factor = 1;
img_processed = zeros(size(img_series));
for frame_iter = 1:numberOfPages
    temp_img = double(img_series(:, :, frame_iter));
    img_processed(:, :, frame_iter) = imgaussfilt(temp_img, gaussian_factor);
end
nucleus_mask = NucleusMask2(img_processed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% read coordinates %%%%%%%%%%%%%%%%%%%%%%%%%
% TrackMate locs
locus_tracks = readmatrix([filepath, filename(1:(end-4)), '.txt'])+1;

resize_factor = 10;
save([filepath, filename(1:(end-4)), '.mat'], "img_series_max", "locus_tracks", "nucleus_mask", "resize_factor", "channel_labels", "pixelSize");
%%
foci_result = struct();
%% dealing with foci channel

for c_iter = 1:2 % DNA and RNA channel

disp(['Processing ', channel_labels{c_iter}, ' channel ...']);

img_series = img_series_max(:, :, :, c_iter);
rc_index = zeros(numberOfPages, 2);

h = size(img_series, 1); w = size(img_series, 2);
numberOfPages = size(img_series, 3);

bw = logical(zeros(size(nucleus_mask)));
if c_iter<=size(locus_tracks, 1)
    bw(locus_tracks(c_iter, 1), locus_tracks(c_iter, 2)) = 1;
    bw = imdilate(bw, true(5, 5));
end

for frame_iter = 1:numberOfPages
    temp_img = imresize(double(imfilter(img_series(:, :, frame_iter), fspecial('average', [3 3]), 'replicate')), resize_factor, "bilinear");
    temp_bw = imresize(bw(:, :, frame_iter), resize_factor, "nearest");
    idx = find(temp_bw);
    if ~isempty(idx)
        [~, i] = max(temp_img(idx));
        [r, c] = ind2sub([h, w]*resize_factor, idx(i));
        rc_index(frame_iter, :) = [r, c]/resize_factor;
    end
end

if rc_index(1, 1)==0 % no RNA
    rc_index = foci_result(1).rc_index;
end

% calculte foci location and intensity
rc_index_rescale = rc_index*resize_factor;
refined_bw = logical(imresize(zeros(size(bw)), resize_factor));
for frame_iter = 1:numberOfPages
    refined_bw(round(rc_index_rescale(frame_iter, 1)), round(rc_index_rescale(frame_iter, 2)), frame_iter) = 1;
end
% TIFwriter(uint8(imdilate(refined_bw, strel('disk', 5))), [filepath, filename(1:(end-4)), '-', channel_labels{c_iter}, '-Center.tif'], 'lzw');

[rna_bkg,base_bkg] = IntensityCalculation(img_series,nucleus_mask, rc_index, 4, 5);
norm_base_bkg = movmean(base_bkg, 10)-500;
intensity = zeros(size(norm_base_bkg));
for frame_iter = 1:numberOfPages
    intensity(frame_iter) = rna_bkg(frame_iter, 1)/norm_base_bkg(frame_iter)*norm_base_bkg(1);
end

if c_iter==1  %  ROI selection
    roi_window = logical(zeros(size(bw)));
    for frame_iter = 1:numberOfPages
        roi_window(ceil(rc_index(frame_iter, 1)), ceil(rc_index(frame_iter, 2)), frame_iter) = 1;
    end
    roi_width = 20;
    roi_resize = imresize(imdilate(roi_window, true(roi_width, roi_width)), resize_factor, "nearest");
end

% export processed DNA and RNA channel images
img_processed_roi = zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages);
img_center = logical(zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages));
temp_refined_bw = imdilate(refined_bw, strel('disk', 3));
for frame_iter = 1:numberOfPages
    disp(['Processing Frame ', num2str(frame_iter), ' ...']);
    temp_img = imresize(double(imgaussfilt(img_series(:, :, frame_iter), 1)), resize_factor);
    temp_mask = temp_refined_bw(:, :, frame_iter);
    temp_roi_resize = roi_resize(:, :, frame_iter);

    [row1, row2, col1, col2] = getROIboundary(temp_roi_resize, roi_width*resize_factor);
    img_processed_roi(row1:row2, col1:col2, frame_iter) = reshape(temp_img(temp_roi_resize), [row2-row1+1, col2-col1+1]);
    img_center(row1:row2, col1:col2, frame_iter) = reshape(temp_mask(temp_roi_resize), [row2-row1+1, col2-col1+1]);
end
TIFwriter(uint16(img_processed_roi), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '.tif']);
TIFwriter(uint8(img_center), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-Center-roi.tif'], 'lzw');

% important parameter: rc_index, bw, intensity
foci_result(c_iter).name = channel_labels{c_iter};
foci_result(c_iter).bw = bw;
foci_result(c_iter).rc_index = rc_index;
foci_result(c_iter).rna_bkg = rna_bkg;
foci_result(c_iter).base_bkg = base_bkg;
foci_result(c_iter).intensity = intensity;
end

save([filepath, filename(1:(end-4)), '.mat'], "foci_result", "roi_window", '-append');

%% dealing with condensate channel
condensate_result = struct; 
%%

for c_iter = 3:4 % OCT4 and BRD4 channel

img_series = img_series_max(:, :, :, c_iter);

% preprocessing parameter

rolling_ball_radius = 50; gaussian_factor = 1;

img_processed = zeros(size(img_series));
for frame_iter = 1:numberOfPages
    temp_img = double(img_series(:, :, frame_iter));
    temp_img = imtophat(temp_img, strel('disk',rolling_ball_radius));
    img_processed(:, :, frame_iter) = imgaussfilt(temp_img, gaussian_factor);
end

roi_width = 20;
roi_bw = imdilate(roi_window, true(roi_width, roi_width))&nucleus_mask;

% HMRF segmentation to define condensate threshold
nclust = 7; seg_point = seg_points(c_iter-2);
HMRFseg_cutoff = zeros(numberOfPages, nclust+1);
for frame_iter = 1:numberOfPages
    temp_img = imresize(img_processed(:, :, frame_iter), resize_factor, "bilinear");
    temp_bw = imresize(roi_bw(:, :, frame_iter), resize_factor, "nearest");
    [HMRFseg, ~] = HMRFseg4img(temp_img, temp_bw, nclust, 0.1, 10^(-8));
    
    for clust_iter = 1:nclust
        HMRFseg_cutoff(frame_iter, clust_iter) = min(HMRFseg.img(HMRFseg.img_class==clust_iter));
    end
    HMRFseg_cutoff(frame_iter, nclust+1) = max(HMRFseg.img(HMRFseg.img_class==nclust));
end

cutoff_new = movmean(HMRFseg_cutoff(:, seg_point), 15);

temp_roi_bw = imdilate(roi_window, true(roi_width, roi_width));
[row1, row2, col1, col2] = getROIboundary(temp_roi_bw, roi_width);
fig1 = figure;
fig1.Units = "inches";
fig1.Position = [7.4,3.9,9.8,6.4];
subplot(1, 2, 1);
imagesc(reshape(img_processed(temp_roi_bw), [row2-row1+1, col2-col1+1]));
colormap(gca, "gray")
daspect([1, 1, 1]);
axis off
subplot(1, 2, 2);
cmap = uint8([4, 0, 0; 56, 46, 142; 137, 48, 141; 215, 31, 40; 239, 127, 25; 244, 191, 27; 244, 237, 70; 255, 255, 255]);
nucleus_partition = reshape(HMRFseg.img_class(imresize(temp_roi_bw, resize_factor, "nearest")), [row2-row1+1, col2-col1+1]*resize_factor);
imagesc(nucleus_partition);
if min(min(nucleus_partition)) == 1
    colormap(gca, cmap(2:end, :));
else
    colormap(gca, cmap);
end
axis off
daspect([1, 1, 1]);
print(fig1, [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-HMRFseg.png'], '-dpng');
close;

% export img-roi
img_processed_roi3 = zeros(roi_width, roi_width, numberOfPages);
for frame_iter = 1:numberOfPages
    temp_roi_bw = imdilate(roi_window(:, :, frame_iter), true(roi_width, roi_width));
    [row1, row2, col1, col2] = getROIboundary(temp_roi_bw, roi_width);
    temp_img_processed = img_processed(:, :, frame_iter);
    img_processed_roi3(:, :, frame_iter) = reshape(temp_img_processed(temp_roi_bw), [row2-row1+1, col2-col1+1]);
end
TIFwriter(uint16(img_processed_roi3), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-processed3.tif']);

% get CD center, interface and boundary, then calculate the distance to
% boundary within roi mask
% roi_width = 20;
roi_resize = imresize(imdilate(roi_window, true(roi_width, roi_width)), resize_factor, "nearest");

CDcenter = logical(zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages));
CDinterface = CDcenter; CDboundary = CDcenter;
img_processed_roi = zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages);
img_processed_roi2 = zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages);
labels = uint16(zeros(roi_width*resize_factor, roi_width*resize_factor, numberOfPages));
for frame_iter = 1:numberOfPages
    disp(['Processing Frame ', num2str(frame_iter), ' ...']);
    temp_img = imresize(img_processed(:, :, frame_iter), resize_factor, "bilinear");
    condensate_bw = temp_img>cutoff_new(frame_iter);
    temp_roi_resize = roi_resize(:, :, frame_iter);
    [row1, row2, col1, col2] = getROIboundary(temp_roi_resize, roi_width*resize_factor);

    temp_CDcenter = getCDcenter(img_processed(:, :, frame_iter), condensate_bw, resize_factor, 0);
    CDcenter(row1:row2, col1:col2, frame_iter) = reshape(temp_CDcenter(temp_roi_resize), [row2-row1+1, col2-col1+1]);
    temp_CDinterface = getCDinterface(temp_img, condensate_bw, 0);
    CDinterface(row1:row2, col1:col2, frame_iter) = reshape(temp_CDinterface(temp_roi_resize), [row2-row1+1, col2-col1+1]);

    img_processed_roi(row1:row2, col1:col2, frame_iter) = reshape(temp_img(temp_roi_resize), [row2-row1+1, col2-col1+1]);
    temp_condensate_bw = logical(zeros(roi_width*resize_factor, roi_width*resize_factor));
    temp_condensate_bw(row1:row2, col1:col2, frame_iter) = reshape(condensate_bw(temp_roi_resize), [row2-row1+1, col2-col1+1]);
    [temp_CDboundary, temp_labels] = getCDboundary(img_processed_roi(:, :, frame_iter), temp_condensate_bw, CDcenter(:, :, frame_iter), CDinterface(:, :, frame_iter), 0);
    CDboundary(:, :, frame_iter) = temp_CDboundary;
    labels(:, :, frame_iter) = temp_labels;

    temp_img = imresize(img_processed(:, :, frame_iter), resize_factor, "nearest");
    img_processed_roi2(row1:row2, col1:col2, frame_iter) = reshape(temp_img(temp_roi_resize), [row2-row1+1, col2-col1+1]);
end

TIFwriter(uint16(img_processed_roi2), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-processed.tif']);
TIFwriter(uint16(img_processed_roi), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-bilinear.tif']);
TIFwriter(uint8(CDcenter), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-CDcenter.tif'], 'lzw');
TIFwriter(uint8(CDinterface), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-CDinterface.tif'], 'lzw');
TIFwriter(uint8(CDboundary), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-CDboundary.tif'], 'lzw');

mask = logical(zeros(size(CDboundary)));
for frame_iter = 1:numberOfPages
    mask(:, :, frame_iter) = img_processed_roi(:, :, frame_iter)>cutoff_new(frame_iter);
end

%%%%%%%%%%%%% export condensate mask %%%%%%%%%%%%%
CDmask = uint8(mask)*50; CDmask(CDboundary)=255;
TIFwriter(uint8(CDmask), [filepath, filename(1:(end-4)), filesep, filename(1:(end-4)), '-', channel_labels{c_iter}, '-CDmask.tif']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dist2bound = zeros(2, numberOfPages);
for f_iter = 1:2 % DNA and RNA channel
    rc_index = foci_result(f_iter).rc_index*resize_factor;
    for frame_iter = 1:numberOfPages
        [row, col] = find(roi_resize(:, :, frame_iter));
        min_row = min(row); min_col = min(col);
        rc_index(frame_iter, :) = rc_index(frame_iter, :)-[min_row-1, min_col-1];
    end
    dist2bound(f_iter, :) = getDist2Bound(rc_index, mask, CDboundary)*pixelSize/resize_factor;
end

condensate_result(c_iter).name = channel_labels{c_iter};
condensate_result(c_iter).img_processed = img_processed;
condensate_result(c_iter).img_processed_roi = img_processed_roi;
condensate_result(c_iter).HMRFseg_cutoff = HMRFseg_cutoff;
condensate_result(c_iter).cutoff = cutoff_new;
condensate_result(c_iter).CDcenter = CDcenter;
condensate_result(c_iter).CDinterface = CDinterface;
condensate_result(c_iter).CDboundary = CDboundary;
condensate_result(c_iter).dist2bound = dist2bound;
condensate_result(c_iter).labels = labels;

end
dist_dna2rna = sqrt(sum((foci_result(1).rc_index-foci_result(2).rc_index).^2, 2))*pixelSize;

% important parameter: condensate_result, dist_dna2rna
save([filepath, filename(1:(end-4)), '.mat'], "condensate_result", "dist_dna2rna", '-append');

%%
close all;
end
