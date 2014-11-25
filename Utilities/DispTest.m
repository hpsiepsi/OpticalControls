function str = uScoreLaTex(str_in)

% find and replace '_' with '\_' so the underscores are correctly displayed
% with a Latex interpreter
 str = strrep(str_in, '_', '\_');


end