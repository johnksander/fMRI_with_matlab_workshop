function [treated_runs,run_TRdata] = treat_runs(data_matrix,TRdata)

%Treat the fmri data run-wise: account for the hemodynamic response,
%linear drift during runs, and standardize voxel timecourses 

run_labels = TRdata.chunks; %could also index as TRdata(:,'chunks') or  TRdata(:,2) 
scan_runs = unique(run_labels); %get a vector of run IDs
num_runs = numel(scan_runs);

treated_runs = cell(num_runs,1); %preallocate a cell aray for the treated data 
run_TRdata = cell(num_runs,1); %and for the treated TR info 


for idx = 1:num_runs
    
    curr_run = scan_runs(idx);
    curr_run = run_labels == curr_run; %  == means "is equal to" and returns a logical 
    run_data = data_matrix(curr_run,:); %use this to index our data matrix 
    run_info = TRdata(curr_run,:); %corresponding TR info 
    
    %---first account for the hemodynamic response 
    %This is typically estimated with a linear model, but we'll use a
    %quick & easy solution that works well for this block-design. Since the
    %TR is 2.5s here, we'll just lag all the data by 2 TRs (roughly 
    %aligning with the peak HDR).
    
    run_data(1:2,:) = []; %assignment with brackets deletes entries 
    run_info = run_info(1:end-2,:); %remove the last two TR entries 
    
    %remove linear drift during the run 
    run_data = detrend(run_data); %super easy 
    
    %zscore each voxel timecourse. This is referred to as "cocktail-blank 
    %normalization" or "mean pattern subtraction", and is the common
    %practice. but see www.frontiersin.org/articles/10.3389/fnins.2013.00174/full
    %for a convicing argument against this! 
    
    run_data = zscore(run_data); 
    
    %put the treated data & info in their cell arrays 
    treated_runs{idx} = run_data; %index cell arrays with curly brackets... mostly 
    run_TRdata{idx} = run_info; 
end

%we really want this data as a matrix, so convert the cell-array 
treated_runs = cell2mat(treated_runs); %we can do this with consistent dims 
run_TRdata = cat(1,run_TRdata{:}); %another way you can do this 









