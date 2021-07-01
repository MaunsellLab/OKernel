function makeFigures(specifierStr)

  if nargin < 1
    specifierStr = 'oneOff';
  end
%   preProcessAll();
  switch specifierStr
    case 'oneOff'
      [T, ~, limits] = getSessionTable('oneOff');
      U = T(T.date >= "2021-05-14" & T.date <= "2021-07-01", :);
      stimProfiles = getOptoProfiles(U);
      plotKernelPage(U, limits, stimProfiles);
    otherwise
      fprintf('makeFigures: unrecognized specifier ''%s''', specifierStr);
  end
end