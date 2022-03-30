function edf1 = edf2mat(top_dir,file_dir,out_dir)

% Inputs:
% top_dir: where you store all your data, e.g., 'Users/Ashley/data/'
% file_dir: the subject/data folder under top_dir, e.g.,'1000/'
% out_dir: the output folder under file_dir, e.g., 'result/'

% output_dir: location to save .mat file

% Outputs:
% eye: structure with all eye tracking data

fprintf('converting .edf file to .mat file...');
fprintf('\n')
tStart = tic;

% get the file name & path
file = dir([top_dir,file_dir,'*.edf']); 

% set up the output directory
out_dir1 = [top_dir,file_dir,out_dir];
if ~isfolder(out_dir1)
    mkdir(out_dir1)
end

for ii = 1:length(file)
% full file path
 full_file_path = [file(ii).folder,'/',file(ii).name];

% check if the file exists
if ~exist(full_file_path)
    error(['Cannot find the file at ', full_file_path])
end

cd(file(ii).folder)
% read in edf file
edf1 = Edf2Mat(file(ii).name);

% save the data file
fprintf('\n');
fprintf('saving .mat file...');
filename = erase(file(ii).name,'.edf');
save([out_dir1,filename,'.mat'],'edf1');
end

tEnd = toc(tStart);
fprintf('\n')
fprintf('done: %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60))
end