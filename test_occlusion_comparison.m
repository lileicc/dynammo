% batch test 

% newton's method

filelist = {'data/c3d/127_07.csv'; ... %run
  'data/c3d/132_32.csv'; ... %bouncy walk
  'data/c3d/132_43.csv'; ... %walk swing shoulder
  'data/c3d/132_46.csv'; ... %walk slow
  'data/c3d/135_07.csv'; ... %mawashigeri
  'data/c3d/141_04.csv'; ... %jump distance
  'data/c3d/141_16.csv'; ... %wave hello
  'data/c3d/143_22.csv'; ... %football throw
  'data/c3d/80_10.csv'};    %boxing

data_names = {'#127_07'; ...
  '#132_32'; ...
  '#132_43'; ...
  '#132_46'; ...
  '#135_07'; ...
  '#141_04'; ...
  '#141_16'; ...
  '#143_22'; ...
  '#80_10'
  };

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
    save(sprintf('temp_result_occlusion_comparison_%d.mat', fileid));
    
    % using Newton direct method
    [X, Y, model, W, bone, bone_var, time, mse] = run_c3d_occlusion_bone(filename, missing_bone, missing_frame_start, missing_frame_end, 'MaxIter', 200, 'Fast', 'NewtonDirect');
    stat{rid} = [stat{rid}; time, mse, missing_frame_start, missing_frame_end, missing_bone];
    save(sprintf('temp_result_occlusion_comparison_%d.mat', fileid));
    
    [X, Y, model, W, bone, bone_var, time, mse] = run_c3d_occlusion_bone(filename, missing_bone, missing_frame_start, missing_frame_end, 'MaxIter', 200, 'Fast', 'SoftBone');
    stat{rid} = [stat{rid}; time, mse, missing_frame_start, missing_frame_end, missing_bone];
    save(sprintf('temp_result_occlusion_comparison_%d.mat', fileid));
    
    % using Newton method
    %[X, Y, model, W, bone, bone_var, time, mse] = run_c3d_occlusion_bone(filename, missing_bone, missing_frame_start, missing_frame_end, 'MaxIter', 200, 'Fast', 'Newton');
    %stat{rid} = [stat{rid}; time, mse, missing_frame_start, missing_frame_end, missing_bone];
    
    % plain bone constraints
    %[X, Y, model, W, bone, bone_var, time, mse] = run_c3d_occlusion_bone(filename, missing_bone, missing_frame_start, missing_frame_end, 'MaxIter', 200, 'Fast');
    %stat{rid} = [stat{rid}; time, mse, missing_frame_start, missing_frame_end, missing_bone];
    
    %tempResult = cell2mat(stat);
    %csvwrite(sprintf('temp_result_occlusion_comparison_%d.csv', fileid), tempResult);
    close all;
  end
  stats{fileid} = stat;  
end


errors = cell(length(filelist), 1);
errors_mean = zeros(length(filelist), 3);
errors_low = errors_mean;
errors_high = errors_mean;
errors_median = errors_mean;
for id = 1 : length(filelist)
  temp = cell2mat(stats{id});
  tmp = reshape(temp(:, 2), 3, []);
  errors{id} = tmp;
  errors_mean(id, :) = mean(tmp, 2)';
  errors_low(id, :) = quantile(tmp', 0.25);
  errors_high(id, :) = quantile(tmp', 0.75);
  errors_median(id, :) = median(tmp, 2)';
end

% plot mean
figure;
colormap(colorGray);
bar(errors_mean);
set(gca, 'XTickLabel', data_names);
legend('dynammo', 'hard constraints', 'soft constraints');
ylabel('average mse');

% boxplot 
% plot error bars
figure;
colormap(colorGray);
bar(errors_median);
set(gca, 'XTickLabel', data_names);
legend('dynammo', 'hard constraints', 'soft constraints');
ylabel('median mse');

hold on;
bar(errors_high, 0.2);
bar(errors_low, 0.5);


figure;
all_errors = cell2mat(errors)';
boxplot(all_errors);


fid = fopen('temp.csv', 'wt');
for id = 1 : length(filelist)
  fprintf(fid, '\n%s\n', filelist{id});
  for rid = 1:REPEATS
    for sid = 1:3
      fprintf(fid, '\t%g\t%g\t%d\t%d\t%d\n', stats{id}{rid}(sid, 1), stats{id}{rid}(sid, 2), stats{id}{rid}(sid, 3), stats{id}{rid}(sid, 4), stats{id}{rid}(sid, 5));
    end
    fprintf(fid, '\n');
  end
end
fclose(fid);