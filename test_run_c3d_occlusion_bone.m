run_c3d_occlusion_bone('data/c3d/132_32.csv', [31, 33], 25, 90, 'MaxIter', 100);

run_c3d_occlusion_bone('data/c3d/132_43.csv', [25, 28], 100, 500, 'MaxIter', 100, 'Fast');


run_c3d_occlusion_bone('data/c3d/80_10.csv', [25, 28, 38, 39, 27], 100, 300, 'MaxIter', 500, 'Fast', 'NewtonDirect');
run_c3d_occlusion_bone('data/c3d/80_10.csv', [25, 28, 38, 39, 27], 100, 300, 'MaxIter', 500, 'Fast', 'Newton');
run_c3d_occlusion_bone('data/c3d/80_10.csv', [25, 28, 38, 39, 27], 100, 300, 'MaxIter', 500, 'Fast');
run_c3d_occlusion('data/c3d/80_10.csv', [25, 28, 38, 39, 27], 100, 300, 'MaxIter', 500, 'Fast');

run_c3d_occlusion_bone('data/c3d/80_10.201-600.csv', [25], 100, 300, 'MaxIter', 1000, 'Fast');


%% 
run_c3d_occlusion_bone('data/c3d/132_43.csv', [25, 28], 100, 500, 'MaxIter', 500, 'Fast', 'NewtonDirect', 'Threshold', 0.04, 'PlotFun', @(x) plot(x([73:75 82:84], :)'));
run_c3d_occlusion_bone('data/c3d/132_43.csv', [25, 28], 100, 500, 'MaxIter', 500, 'Fast', 'Newton', 'Threshold', 0.04);
run_c3d_occlusion_bone('data/c3d/132_43.csv', [25, 28], 100, 500, 'MaxIter', 500, 'Fast', 'Threshold', 0.04);
run_c3d_occlusion('data/c3d/132_43.csv', [25, 28], 100, 500, 'MaxIter', 500, 'Fast');


%% occlude all right hand markers
run_c3d_occlusion_bone('data/c3d/132_43.csv', [25, 28, 38, 39, 27], 100, 500, 'MaxIter', 500, 'Fast', 'NewtonDirect');
run_c3d_occlusion_bone('data/c3d/132_43.csv', [25, 28, 38, 39, 27], 100, 500, 'MaxIter', 500, 'Fast', 'NewtonDirect', 'Threshold', 0.04);
run_c3d_occlusion_bone('data/c3d/132_43.csv', [25, 28, 38, 39, 27], 100, 500, 'MaxIter', 500, 'Fast', 'Newton', 'Threshold', 0.04);
run_c3d_occlusion_bone('data/c3d/132_43.csv', [25, 28, 38, 39, 27], 100, 500, 'MaxIter', 500, 'Fast', 'Threshold', 0.04);
run_c3d_occlusion('data/c3d/132_43.csv', [25, 28, 38, 39, 27], 100, 500, 'MaxIter', 500, 'Fast');
