function OPExamples()

  fileNames = {'2019-02-25',  '2019-02-15', '2019-02-05'};
  addpath('/Users/Shared/Data/OKernel/ Analysis/Analysis Code');
  fH = figure(1);
  set(fH, 'units', 'inches', 'position', [26.5, 7, 8.5, 11]);
  clf;
  for f = 1:length(fileNames)
    load(sprintf('645/Matfiles/%s.mat', fileNames{f}), 'file', 'trials');
    OPAMatlab(file, trials, f);
  end
end  