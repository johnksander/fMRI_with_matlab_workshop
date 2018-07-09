clear          %remove all variables from the workspace
clc            %clears the command window
format compact %removes extra spaces from the command window's output 
%always start scripts with at least a clear!


%----specify the directories needed in your analysis----
%I always use absolute directories, this avoids unintended consequences with
%duplicate files, mixing up working directories, etc. 
%I often start by hardcoding one "base" directory, then building my
%directory tree out from there. 

base_dir = '~/Desktop/work/fMRI_with_matlab_workshop';
data_dir = fullfile(base_dir,'dataset'); %fullfile() is an easy way to create directories 
func_dir = fullfile(base_dir,'helper_functions');
NTB_dir = fullfile(base_dir,'nifti_toobox');
%If you want to use any extra functions or libraries, they must be in
%matlab's path. Matlab won't know they exist otherwise 
addpath(func_dir) %add our helper_functions directory to the path
addpath(NTB_dir)
%side note on the nifti toolbox: matlab now has native .nii functionality...
%I just havn't tried it out yet www.mathworks.com/help/images/ref/niftiread.html

%----load and inspect a subject's anatomical volume----
    
subj_anat = fullfile(data_dir,'subj_101','anat.nii.gz');
%don't worry about the .gz at the end, this is common & handled by NTB

subj_anat = load_nii(subj_anat);
%Nifti files contain a lot of metadata as well. We'll ignore that for now,
%we just want the image volume (3-D scan matrix). This file was loaded as a
%structure datatype, with subfields for scans & metadata. 

subj_anat = subj_anat.img; %index the "img" subfield for the scan volume
%we now have the anatomical volume as a 3-D matrix (124 x 256 x 256)
%it's good practice to note your scan dimensions, we'll revisit this later 

%let's visualize an axial slice 
anat_slice = subj_anat(:,:,135); %index matricies with parentheses 
imagesc(anat_slice) %visualize 2-D slice, check for brain. 

%----loop through subjects & preprocess functional data----
subjects = 101:103; %create a vector of our 3 subjects' ID numbers
num_subs = numel(subjects);
subj_preproc_data = cell(num_subs,1); %initialize cell array for the preprocessed data  
subj_class_labels = cell(num_subs,1); %and for the class & run labels while we're at it 
subj_run_labels = cell(num_subs,1); 


for idx = 1:num_subs %loop over subjects
    
    curr_subject = subjects(idx); %current subject ID #
    %specifiy the subject's data directory
    subj_dir = ['subj_' num2str(curr_subject)]; %brackets concatenate like-datatypes
    subj_dir = fullfile(data_dir,subj_dir); 
    
    %load the subject's functional data
    fmri_data = load_nii(fullfile(subj_dir,'bold.nii.gz')); %this will take a moment
    fmri_data = fmri_data.img;
    
    size(fmri_data) %look at the matrix dimensions (40 x 64 x 64 x 1452) 
    %that's a 4-D matrix containing 1452 scans (each size 40 x 64 x 64)
    
    %load the subject's ROI mask 
    ROI_mask = load_nii(fullfile(subj_dir,'mask4_vt.nii.gz'));
    ROI_mask = ROI_mask.img; 
    ROI_mask = logical(ROI_mask); %convert to logical (boolean) datatype 
    size(ROI_mask) %ensure this volume's size matches your data 
    
    %now we need to extract the ROI data, and flatten the 4-D timeseries
    %into a 2-D matrix for analysis where rows are timepoints (scans) 
    %and columns are voxels. Try making this function!
    
    data_matrix = mask_data(fmri_data,ROI_mask);
    
    %load the TR file 
    %side note: there's many ways to deal with reading text files in matlab
    %heres a helpful guide www.mathworks.com/help/matlab/import_export/ways-to-import-text-files.html

    TRdata = readtable(fullfile(subj_dir,'labels.txt'));
    %this gives us a table datatype, with a row for each scan. Columns
    %contain the stimulus labels, and run numbers (called "chunks" here?). 
    
    [data_matrix,TRdata] = treat_runs(data_matrix,TRdata);
    
    subj_preproc_data{idx} = data_matrix;
    subj_class_labels{idx} = TRdata.labels;
    subj_run_labels{idx} = TRdata.chunks;
    
    %add a little message so we can see progress 
    fprintf('finished preprocessing subject %i \n',curr_subject) 
end


%now lets set up the decoding analysis, let's do this within-subject as well 

subj_accs = NaN(num_subs,1); %for our decoding accuracies 

for idx = 1:num_subs %loop over subjects
    
    data_matrix = subj_preproc_data{idx};
    class_labels = subj_class_labels{idx};
    run_labels = subj_run_labels{idx};
    
    %we're not interested in the "resting" trials, so let's get rid of them
    rest_trials = ismember(class_labels,'rest'); %gives us a logical 
    
    %use the tilde to indicate "is not" 
    data_matrix = data_matrix(~rest_trials,:); 
    class_labels = class_labels(~rest_trials);
    run_labels = run_labels(~rest_trials);
    
    
    %let's do a leave-one-run-out cross-validation scheme 
    %In all decoding analyses, you must maintain strict independence between 
    %the training and testing sets. In the context of fMRI, this means data
    %from a given run must be entirely commited to either the testing or
    %training set, but never both (e.g. if the triaining set contains data
    %from run #3, the testing set cannot include any data from run #3).
    %This is also why preprocessing is perfomed within-run, so there are no
    %dependencies in the data aross runs. 
    
    
    scan_runs = unique(run_labels); %get a vector of run IDs
    num_runs = numel(scan_runs);
    CV_accs = NaN(num_runs,1); %for our CV accuracies 
    
    for run_idx = 1:num_runs %loop over runs (cross-validation loop) 
        
        testing_run = scan_runs(run_idx); 
        %get logicals for individual trials 
        testing_trials = run_labels == testing_run;
        training_trials  = ~testing_trials; 
        
        
        %divy up the data
        training_data = data_matrix(training_trials,:); 
        training_lables = class_labels(training_trials);
        %indexing a cell-array with parenthases returns cells, rather than
        %their contents (as with {} indexing) 
        
        testing_data = data_matrix(testing_trials,:);
        testing_labels = class_labels(testing_trials,:);
        
        %train a linear discriminant analysis model 
        fit_mdl = fitcdiscr(training_data,training_lables);
        %test the model 
        predictions = predict(fit_mdl,testing_data);
        %check and see how accurate
        correct_preds = strcmpi(testing_labels,predictions);
        %record average accuracy
        CV_accs(run_idx) = sum(correct_preds) / numel(correct_preds);
    
    end
    
    subj_accs(idx) = mean(CV_accs); %average accuracy across CV folds
    
end

disp(subj_accs)

%if we have extra time, we can do one of these: 
%1) evaluate the resutls with permutation testing 
%2) visualize the fMRI data with PCA
%3) HDR estimation 




