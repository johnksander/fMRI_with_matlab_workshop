function data_matrix = mask_data(bold_data,mask)
%inputs: 4-D bold timeseries, 3-D anatomical mask
%output: 2-D analysis matrix of ROI data


num_voxels = sum(mask(:)); %number of voxels in the ROI
num_scans = size(bold_data,4); %number of scans (timepoints) in the bold series

%preallocate the flattened matrix (timepoints x voxels)
data_matrix = NaN(num_scans,num_voxels); 

for scanidx = 1:num_scans %loop over scans
    curr_scan = bold_data(:,:,:,scanidx);
    curr_scan = curr_scan(mask); %select ROI voxels with logical matrix 
    data_matrix(scanidx,:) = curr_scan; %assign the voxel vector to a row
end

