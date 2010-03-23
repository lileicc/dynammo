% batch test 

% newton's method

filelist = {'data/c3d/127_07.csv'; ...
'data/c3d/132_32.csv'; ...
'data/c3d/132_43.csv'; ...
'data/c3d/132_46.csv'; ...
'data/c3d/135_07.csv'; ...
'data/c3d/141_04.csv'; ...
'data/c3d/141_16.csv'; ...
'data/c3d/143_22.csv'; ...
'data/c3d/80_10.csv'};

REPEATS = 10;
stats = cell(1, length(filelist));
stat = cell(REPEATS, 1);

for fileid = 1:length(filelist)
  
  filename = filelist{fileid};
  lab = importdata(filename);
  N = size(lab.data, 1);
  M = size(lab.data, 2); 
  for rid = 1 : REPEATS
    stat{rid} = [];   
    
    missing_length = unidrnd(N);
    missing_frame_start = unidrnd(N - missing_length + 1);
    missing_frame_end = missing_frame_start + missing_length - 1;
    missing_bone = unidrnd(round(M/3));
    
    % using DynaMMo
    [X, Y, model, W, bone, bone_var, time, mse] = run_c3d_occlusion(filename, missing_bone, missing_frame_start, missing_frame_end, 'MaxIter', 200, 'Fast');
    stat{rid} = [stat{rid}; time, mse, missing_frame_start, missing_frame_end, missing_bone];
    
    % using Newton direct method
    [X, Y, model, W, bone, bone_var, time, mse] = run_c3d_occlusion_bone(filename, missing_bone, missing_frame_start, missing_frame_end, 'MaxIter', 200, 'Fast', 'NewtonDirect');
    stat{rid} = [stat{rid}; time, mse, missing_frame_start, missing_frame_end, missing_bone];
    
    % using Newton method
    [X, Y, model, W, bone, bone_var, time, mse] = run_c3d_occlusion_bone(filename, missing_bone, missing_frame_start, missing_frame_end, 'MaxIter', 200, 'Fast', 'Newton');
    stat{rid} = [stat{rid}; time, mse, missing_frame_start, missing_frame_end, missing_bone];
    
    % plain bone constraints
    [X, Y, model, W, bone, bone_var, time, mse] = run_c3d_occlusion_bone(filename, missing_bone, missing_frame_start, missing_frame_end, 'MaxIter', 200, 'Fast');
    stat{rid} = [stat{rid}; time, mse, missing_frame_start, missing_frame_end, missing_bone];
    
    tempResult = cell2mat(stat);
    csvwrite(sprintf('temp_result_occlusion_comparison_%d.csv', fileid), tempResult);
    close all;
  end
  stats{fileid} = stat;
  
  
end
