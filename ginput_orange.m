function varargout = ginput_orange(arg1)
%GINPUT_ORANGE Graphical input from mouse.
%   [X,Y] = GINPUT_ORANGE(N) gets N points from the current axes and returns
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse.  Data points are entered by pressing
%   a mouse button or any key on the keyboard except carriage return,
%   which terminates the input before N points are entered.
%
%   [X,Y] = GINPUT_ORANGE gathers an unlimited number of points until the
%   return key is pressed.
%
%   [X,Y,BUTTON] = GINPUT_ORANGE(N) returns a third result, BUTTON, that
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
%
%   Examples:
%       [x,y] = ginput_orange;
%
%       [x,y] = ginput_orange(5);
%
%       [x, y, button] = ginput_orange(1);
%
%   See also GTEXT, WAITFORBUTTONPRESS.

%   Copyright 1984-2013 The MathWorks, Inc.

x = []; y = []; userInputs = [];

if ~matlab.ui.internal.isFigureShowEnabled
    error(message('MATLAB:hg:NoDisplayNoFigureSupport', 'ginput_orange'))
end
    
    % Check inputs
    if nargin == 0
        how_many = -1;
    else
        how_many = arg1;
        
        if ~isPositiveScalarIntegerNumber(how_many)
            error(message('MATLAB:ginput_orange:NeedPositiveInt'))
        end
        if how_many == 0
            % If input argument is equal to zero points,
            % give a warning and return empty for the outputs.
            warning (message('MATLAB:ginput_orange:InputArgumentZero'));
        end
    end
    
    % Get figure
    fig = gcf;
    figure(fig);
    
    % Make sure the figure has an axes
    gca(fig);
    
    % Setup the figure to disable interactive modes and activate pointers. 
    initialState = setupFcn(fig);
    
    % onCleanup object to restore everything to original state in event of
    % completion, closing of figure errors or ctrl+c. 
    c = onCleanup(@() restoreFcn(initialState));
       
    drawnow
    char = 0;
    
    while how_many ~= 0
        try
            mode = waitForUserInput(fig);
        catch %#ok<CTCH>
            cleanup(c);
            if(ishghandle(fig))
                error(message('MATLAB:ginput_orange:Interrupted'));
            else
                error(message('MATLAB:ginput_orange:FigureDeletionPause'));
            end
        end
        
        
        % Make sure figure has not been closed
        checkFigureAvailable();

        if (isCorrectFigure(fig))
            switch mode
                case 'key'
                    char = get(fig, 'CurrentCharacter');
                    curUserInput = abs(get(fig, 'CurrentCharacter'));
                case 'mouse'
                    curUserInput = get(fig, 'SelectionType');
                    if strcmp(curUserInput,'open')
                        curUserInput = 1;
                    elseif strcmp(curUserInput,'normal')
                        curUserInput = 1;
                    elseif strcmp(curUserInput,'extend')
                        curUserInput = 2;
                    elseif strcmp(curUserInput,'alt')
                        curUserInput = 3;
                    else
                        error(message('MATLAB:ginput_orange:InvalidSelection'))
                    end
            end
            axes_handle = get(fig,'CurrentAxes');
            pt = get(axes_handle, 'CurrentPoint');

            how_many = how_many - 1;

            if(char == 13)
                % if the return key was pressed, char will == 13,
                % and that's our signal to break out of here whether
                % or not we have collected all the requested data
                % points.
                % If this was an early breakout, don't include
                % the <Return> key info in the return arrays.
                % We will no longer count it if it's the last input.
                break;
            end

            x = [x;pt(1,1)]; %#ok<AGROW>
            y = [y;pt(1,2)]; %#ok<AGROW>
            userInputs = [userInputs;curUserInput]; %#ok<AGROW>
        end
    end
    
    % Cleanup and Restore 
    cleanup(c);
    
    if nargout == 1
        varargout{1} = [x y];
    else
        varargout{1} = x;
    end
    if nargout > 1
        varargout{2} = y;
    end
    if nargout > 2
        varargout{3} = userInputs;
    end

    
end

function valid = isPositiveScalarIntegerNumber(how_many)
 
valid = ~ischar(how_many) && ...            % is numeric
        isscalar(how_many) && ...           % is scalar
        (fix(how_many) == how_many) && ...  % is integer in value
        how_many >= 0;                      % is positive
end

function mode = waitForUserInput(fig)
waitfor(fig,'UserData')
% Extract mode to determine if key or mouse was used
mode = get(fig,'UserData'); 
if ischar(mode)
    ud = strsplit(mode, '_');
    mode = ud{1};
end% Reset user data to prepare for next trigger
set(fig,'UserData',[])
end

function initialState = setupFcn(fig)

% Store Figure Handle. 
initialState.figureHandle = fig; 

% Suspend figure functions
initialState.uisuspendState = uisuspend(fig);

% Disabling ^C for edit menu so the only ^C is for interrupting the function
initialState.AcceleratorMenu = findall(fig,'Type','uimenu','Accelerator','C');
set(initialState.AcceleratorMenu,'Accelerator','');

% Extract user data
initialState.PreviousUserData = get(fig,'UserData');

% Set callbacks to distinguish from key and mouse triggers
% Using random numbers to distinguish things like double clicks.
set(fig, 'WindowButtondownFcn', @(~,~)set(fig, 'UserData', ['mouse_' num2str(rand)]));
set(fig, 'WindowKeyPressFcn', @(~,~) set(fig, 'UserData', ['key_' num2str(rand)]));

% Disable Plottools Buttons
initialState.toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn')];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end

%Setup empty pointer
cdata = NaN(16,16);
hotspot = [8,8];
set(gcf,'Pointer','custom','PointerShapeCData',cdata,'PointerShapeHotSpot',hotspot)

% Create uicontrols to simulate fullcrosshair pointer.
initialState.CrossHair = createCrossHair(fig);

% Adding this to enable automatic updating of currentpoint on the figure 
% This function is also used to update the display of the fullcrosshair
% pointer and make them track the currentpoint.
set(fig,'WindowButtonMotionFcn',@(o,e) updateCrossHair(o,initialState.CrossHair));

% Get the initial Figure Units
initialState.fig_units = get(fig,'Units');
end

function restoreFcn(initialState)
if ishghandle(initialState.figureHandle)
    delete(initialState.CrossHair);
    
    % Figure Units
    set(initialState.figureHandle,'Units',initialState.fig_units);
    
    % Reset user data
    set(initialState.figureHandle,'UserData',initialState.PreviousUserData)
    
    % Enable Ctrl+c
    set(initialState.AcceleratorMenu,'Accelerator','C');
    
    % Plottools Icons
    if ~isempty(initialState.toolbar) && ~isempty(initialState.ptButtons)
        set (initialState.ptButtons(1),'Enable',initialState.ptState{1});
        set (initialState.ptButtons(2),'Enable',initialState.ptState{2});
    end
    
    % UISUSPEND
    uirestore(initialState.uisuspendState);    
end
end

function updateCrossHair(fig, crossHair)
% update cross hair for figure.
gap = 3; % 3 pixel view port between the crosshairs
cp = hgconvertunits(fig, [fig.CurrentPoint 0 0], fig.Units, 'pixels', fig);
cp = cp(1:2);
figPos = hgconvertunits(fig, fig.Position, fig.Units, 'pixels', fig.Parent);
figWidth = figPos(3);
figHeight = figPos(4);

% Early return if point is outside the figure
if cp(1) < gap || cp(2) < gap || cp(1)>figWidth-gap || cp(2)>figHeight-gap
    return
end

set(crossHair, 'Visible', 'on');
thickness = 1; % 1 Pixel thin lines. 
set(crossHair(1), 'Position', [0 cp(2) cp(1)-gap thickness]);
set(crossHair(2), 'Position', [cp(1)+gap cp(2) figWidth-cp(1)-gap thickness]);
set(crossHair(3), 'Position', [cp(1) 0 thickness cp(2)-gap]);
set(crossHair(4), 'Position', [cp(1) cp(2)+gap thickness figHeight-cp(2)-gap]);
end

function crossHair = createCrossHair(fig)
% Create thin uicontrols with black backgrounds to simulate fullcrosshair pointer.
% 1: horizontal left, 2: horizontal right, 3: vertical bottom, 4: vertical top
for k = 1:4
    crossHair(k) = uicontrol(fig, 'Style', 'text',...
                             'Visible', 'off',...
                             'Units', 'pixels',...
                             'BackgroundColor', [1 0.5 0],...
                             'HandleVisibility', 'off',...
                             'HitTest', 'off'); %#ok<AGROW>
end
end

function checkFigureAvailable()
% See if root has children
figchildren = allchild(0);
if isempty(figchildren)
    error(message('MATLAB:ginput_orange:FigureUnavailable'));
end
end

function valid = isCorrectFigure(fig)
% g467403 - ginput_orange failed to discern clicks/keypressed on the figure it was 
% registered to if the figure's handleVisibility was set to 'callback'
figchildren = allchild(0);
% Select figure at top
ptr_fig = figchildren(1);
% Check if they match
valid = isequal(ptr_fig,fig);
end

function cleanup(c)
if isvalid(c)
    delete(c);
end
end