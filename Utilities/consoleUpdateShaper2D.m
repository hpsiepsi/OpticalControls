function dispDat = consoleUpdateShaper2D(dispDat, zoom2D, w_e, w_d)
% look for the GUI in the baseworkspace
if (dispDat.gui == 0)
    % no GUI so make new
    propName = {'HorizontalAlignment'};
    propVal  = {'center'};
    boxVal   = {'center'};
    pnlpVal  = {'center'};


    %prog      = ProgressImage(0, 1.0, expProgPhase);
    firstCol  = 10;
    txtWidth  = 100;
    firstBox  = firstCol + txtWidth + 10;
    boxWidth  = 120;
    secondCol = firstCol + txtWidth  + boxWidth + 20;
    secondBox = secondCol + txtWidth + 10;
    
    dispDat.gui = figure(floor(rand(1,1)*1e9));
    
    set(dispDat.gui, 'MenuBar', 'none', 'NumberTitle','off', ...
        'Name', 'Monitor 2DIR Shaper Analysis','Position',[250 250 700 500], ...
        'Toolbar', 'none', 'Resize', 'off');

    
    % make text boxes
    uicontrol('Parent',dispDat.gui,'Style', 'text', 'Position',[firstCol 440 txtWidth 20],'String', ...
        'Processing file ', propName,pnlpVal);
    dispDat.hd_fileNum = uicontrol('Parent',dispDat.gui,'Style', 'edit', 'Position',[firstBox 440 boxWidth 25],'String', ...
        '',propName,boxVal);
    dispDat.hd_fileName = uicontrol('Parent',dispDat.gui,'Style', 'edit', 'Position',[firstCol 410 (boxWidth+txtWidth+10) 25],'String', ...
        '',propName,boxVal);

    % Status
    uicontrol('Parent',dispDat.gui,'Style', 'text', 'Position',[firstCol 470 txtWidth-50 20],'String', ...
        'Status ', propName,pnlpVal);
    dispDat.hd_status = uicontrol('Parent',dispDat.gui,'Style', 'edit', 'Position',[firstBox-50 470 boxWidth+50 25],'String', ...
        '',propName,boxVal);
    
    %Time
    uicontrol('Parent',dispDat.gui, 'Style', 'text', 'Position',[firstCol 380 txtWidth 20],'String', ...
        'Time so far ', propName,pnlpVal);
    dispDat.hd_tTime = uicontrol('Parent',dispDat.gui, 'Style', 'edit', 'Position',[firstBox 380 boxWidth 25],'String', ...
        '0 min', propName,boxVal);
    uicontrol('Parent',dispDat.gui, 'Style', 'text', 'Position',[firstCol 350 txtWidth 20],'String', ...
        'Time last file ', propName,pnlpVal);
    dispDat.hd_lTime = uicontrol('Parent',dispDat.gui, 'Style', 'edit', 'Position',[firstBox 350 boxWidth 25],'String', ...
        '0 sec', propName,boxVal);

    % Noise and amplitude
    uicontrol('Parent',dispDat.gui, 'Style', 'text', 'Position',[firstCol 320 txtWidth 20],'String', ...
        'Noise Last File', propName,pnlpVal);
    dispDat.hd_noise = uicontrol('Parent',dispDat.gui, 'Style', 'edit', 'Position',[firstBox 320 boxWidth 25],'String', ...
        '0 %', propName,boxVal);
    uicontrol('Parent',dispDat.gui, 'Style', 'text', 'Position',[firstCol 290 txtWidth 20],'String', ...
        'Amplitude Last File', propName,pnlpVal);
    dispDat.hd_amp = uicontrol('Parent',dispDat.gui, 'Style', 'edit', 'Position',[firstBox 290 boxWidth 25],'String', ...
        '0', propName,boxVal);
    % Axes
    dispDat.hd_axes = axes('Units','pixels','Parent', dispDat.gui, 'Position', [270, 50, 400, 400]);
    mesh(peaks(20))
end

% a 6 line output is used to display several data items, define the lines
    set(dispDat.hd_fileNum,  'String', [ num2str(dispDat.curNum) ' of ' num2str(dispDat.tNum)]);
    set(dispDat.hd_fileName, 'String', dispDat.filename);
    set(dispDat.hd_tTime,    'String', [num2str(dispDat.tTime/60) ' min'])
    set(dispDat.hd_lTime,    'String', [num2str(dispDat.tLast) ' sec'])  
    set(dispDat.hd_noise,    'String', [num2str(dispDat.noiseLast) ' %']) 
    set(dispDat.hd_amp,      'String', [num2str(dispDat.ampLast, '%10.2e') ' ']) 
    set(dispDat.hd_status,   'String', dispDat.status)
    
    if ~(isempty(zoom2D))
%       contourf(dispDat.hd_axes,w_e, w_d, zoom2D',40, 'LineStyle', 'none'); 
%       line([w_e(1), w_e(end)],[w_e(1), w_e(end)], 'LineStyle', ':', 'Color', [.7 .7 .7], 'LineWidth', 4.0)
    end
    
%    if (dispDat.tNum > 0)&&(dispDat.curNum == dispDat.tNum)
%        set(dispDat.gui,'Name', 'Monitor 2DIR Shaper Analysis (Complete)' )
%    end
    
    drawnow
    
end