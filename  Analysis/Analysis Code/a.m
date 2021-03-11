function a

  dirName = '/Users/Shared/Data/OKernel/902/MatFiles';
  dirStructs = dir(dirName);                              % data directory contents
  fileNames = {dirStructs(:).name};                       % data directory file names
  numFiles = length(fileNames);
  for f = 1:numFiles
      validFiles(f) = isstrprop(fileNames{f}(1), 'digit') && length(fileNames{f}) > 1; %#ok<AGROW>
  end
  fileNames = {fileNames{validFiles}};
%   fileValues = str2double(fileNames);             % sort the files numerically
%   [~, indices] = sort(fileValues);
%   fileNames = {fileNames{indices}};

  for f = 1:numFiles
    load(fullfile(dirName, fileNames{f}), 'file', 'trials');
    OKMatlab([], file, trials);
  end  
end