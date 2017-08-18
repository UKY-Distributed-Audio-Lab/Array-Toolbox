function varargout = ListenBeamform(varargin)
% LISTENBEAMFORM MATLAB code for ListenBeamform.fig
% This code constructs the GUI which allows the user to select and listen to 
% detected streams found during ASA. Users have the ability to select 
% streams and listen to their beamformed wave files.
% Users may select which beamformed files to perform the masking component
% on.  This step is provided to allow for more efficient masking on reliable
% streams, rather than masking both dominant streams and their sidelobes.
%
% This code produces a figure illustrating masked sound sources on a 2D
% plot. The starting point and ending point of each source is shown and
% each are labeled with a number. A separate figure is produced which
% displays two list boxes with the found streams in the list box labeled 
% 'Detected Streams'. The two pushbuttons - play and stop - located beside
% this list box allow users to listen to beamformed wave files associated 
% with the detected streams labeled in the box. The user has the ability to
% add beamformed streams from the left list box to the list box on the
% right which is labeled 'Selected Streams to Mask', as well as remove 
% items in the right box. The user can then select the button 'Perform Masking'
% to perform the masking component only on items added to the right list box.
%
%
% Written by Kirstin Brangers       July 2012
% Last Modified                     15-Aug-2012



% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ListenBeamform_OpeningFcn, ...
                   'gui_OutputFcn',  @ListenBeamform_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% --- Executes just before ListenBeamform is made visible.
function ListenBeamform_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ListenBeamform (see VARARGIN)


% Load data to be used in GUI

% From Mainscript.m:
% scriptout = [Xstart, Ystart, Zstart, Xend, Yend, Zend, StartTime, Filename]
%           scriptout is a cell array
% stseg = [Segment/Stream Number, Intensity, X, Y, Z, Frame Number, Detection Indicator]
% fovp structure => has mic positions


% Load Script with Masked Data       
bffname = varargin{1};              % Filename of Beamformed Data Script Files
bffname = char(bffname);            % Convert cell array to string
load(bffname);                      % Load data
handles.bffname = bffname;          % Save filename in handles

% load CocktailClip3Secs_BothTest_ManualEntrybfscript.mat

handles.basenam = basenam;
handles.pn = []
handles.movieflag = moviemake_flag;
handles.trez = fovp.trez;
handles.scriptout = scriptout;
handles.stseg = stseg;
handles.fovp = fovp;

% Load Start and End Positions (from scriptout)
handles.Xstart = cell2mat(scriptout(:,2));    handles.Ystart = cell2mat(scriptout(:,3));    handles.Zstart = cell2mat(scriptout(:,4));
handles.Xend = cell2mat(scriptout(:,5));      handles.Yend = cell2mat(scriptout(:,6));      handles.Zend = cell2mat(scriptout(:,7));

% Get Screen Size for Plotting Figures
monitorP = get(0,'Screensize');         % Primary Monitor
dual = get(0,'MonitorPositions');
[mon c] = size(dual);
if mon==2
    monitor2 = [dual(2,1), dual(2,2), dual(2,3)-dual(2,1), dual(2,4)];    % Secondary Monitor
    w = min(monitor2(3),monitorP(3));
    h = min(monitor2(4),monitorP(4));
    screen = [1 1 w h];
else
    screen = monitorP;
end

    
% Load Start Times of Segments
handles.StartTimes = cell2mat(scriptout(:,8));                

% Find Number of Streams Detected - Max Number of Frames
[streams cols] = size(scriptout);

% Check that length of X,Y,Z vectors equal number of streams
if (length(handles.Xstart)~=streams || length(handles.Ystart)~=streams || length(handles.Zstart)~=streams)
    error('Number of segments does not equal maximum number of Start streams found')
end
if (length(handles.Xend)~=streams || length(handles.Yend)~=streams || length(handles.Zend)~=streams)
    error('Number of segments does not equal maximum number of End streams found')
end


% Load WAVE files associated with streams after masking and beamforming
% Store data from filename into matrix Y
% Each column is a seperate stream,  # columns = # streams
% WAVE file is stored in each column of Y
% Sample rate is stored in each column of FS
handles.pns = handles.pn;
for i=1:streams
    filename = scriptout{i,9};    % Cycle through filenames
    fname = char(filename);     % Convert from cell array
    fstatus = exist([handles.pns, fname],'file');
    if fstatus == 2
        [y(:,i),fs(:,i)] = wavread([handles.pns, fname]);  % Load signal data
    else 
        [fname, handles.pns] = uigetfile('*.wav',['Find created wavefile for stream ' int2str(i)]);
        [y(:,i),fs(:,i)] = wavread([handles.pns, fname]);  % Load signal data        
    end
end


% Find max number of samples in Y to determine maximum time
handles.maxY=max(length(y(:,:)));
maxFS = max(fs);
handles.maxTime = handles.maxY/maxFS;        % Max length/time of sound clip


% Update number of streams
handles.leftstreams = streams;

% Update Y and FS
handles.y(:,:) = y(:,:);
handles.fs(1,:) = fs(1,:);

% Compute Limits for Plots
minx = min(fovp.mp(1,:));
maxx = max(fovp.mp(1,:));
miny = min(fovp.mp(2,:));
maxy = max(fovp.mp(2,:));
minz = min(fovp.mp(3,:));
minz = min([minz,0]);
maxz = max(fovp.mp(3,:));

% Figure - Detected Sound Sources
posSourceFig = [30, 150, 1000, 800];            % Position of Figure
% Check that figure will fit on screen
if posSourceFig(3)+posSourceFig(1)>screen(3)	% Width               
    posSourceFig(3)=screen(3);
    posSourceFig(1)=1;
    posSourceFig(2)=1;
end
if posSourceFig(4)+posSourceFig(2)>screen(4)	% Height
    posSourceFig(4)=screen(4);
    posSourceFig(1)=1;
    posSourceFig(2)=1;
end


% System Background Color
bkcolor = get(0,'DefaultUicontrolBackgroundColor');
% Source figure
handles.hSources = figure('Name','Detected Streams',...
                'NumberTitle','off',...
                'Color', bkcolor,...
                'Position',posSourceFig);
axis off
hold on
cindx = 1;                          % Color/Marker index
% Plot Detected Sound Sources in 2D Plot
for st=1:streams
    colors = ['r' 'b' 'm' 'k','g'];     % Available source colors
    markers = ['.' '.' '.' '.' '.'];    % Default markers 
    len = length(colors);               % Length of matrix holding available colors
    S = 30;                             % Size of source labels
    
    % Plot points - Start to End
    % Start
    gscatter(handles.Xstart(st,1),handles.Ystart(st,1),st,colors(cindx),markers(cindx),S,'off', 'X Meters', 'Y Meters');   % Starting Points
    text(handles.Xstart(st,1)+.04,handles.Ystart(st,1),[num2str(st),'Start'],'FontSize',8);         % Source Labels - Source#start
    % End
    gscatter(handles.Xend(st,1),handles.Yend(st,1),st,colors(cindx),markers(cindx),S);              % Ending Points
    text(handles.Xend(st,1)+.04,handles.Yend(st,1),[num2str(st)],'FontSize',8);                     % Source Labels - Source#end
    
    % Connect Start and End points with solid black line
    plot([handles.Xstart(st,1) handles.Xend(st,1)],[handles.Ystart(st,1) handles.Yend(st,1)],'-k');
    
    % Properties of Plot
    title('Detected Streams - Start and End Positions')
    grid on
    axis on
    axis([minx-0.5,maxx+0.5,miny-0.5,maxy+0.5])
    legend off
    
    % Increment Counter
    cindx = cindx+1;
    if cindx > len
        cindx = 1;
    end
    
end
pause(1)
hold off

% Create a List Box that has available sources
% Users can select the source and then press play button
% Lengths of source names must be the same to prevent error from occurring
handles.sSources = [];
str = zeros(1,streams);
for st=1:streams                            % ALL beginning available streams
    str(1,st) = st;
    if st<10
        handles.sSources = [handles.sSources; ['Stream0', num2str(st)]];    
    elseif st<100
        handles.sSources = [handles.sSources; ['Stream', num2str(st)]];
    else
        disp('100 streams have been detected. Change secondary thresholds to minimize number of streams detected');
    end        
end


% Set up listbox to include all found sources
set(handles.lb_sourcesLeft, 'String', handles.sSources);

% Create defualt audioplayers - Default: Source1
handles.audio1 = audioplayer(handles.y(:,1),handles.fs(1,1));

% Default stream: Stream1
handles.sources = 1;                % Default stream to play audio for
handles.leftstreams = str;          % Available streams in Left LB    
handles.allstreams = str;           % All streams found
handles.Leftstrings = handles.sSources; % String of streams available in Left LB
% Initialize streams
handles.maskstreams = [];           % Streams in Right LB
handles.Maskstrings = [];           % Strings of streams in Right LB

% Set position/location of figure
posFigure = get(hObject,'Position');                % Get position/size of figure from GUIDE
position = [225, 50, posFigure(3), posFigure(4)];   % Set position of figure
set(hObject,'Position', position);

% Choose default command line output for ListenBeamform
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ListenBeamform wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% --- Outputs from this function are returned to the command line.
function varargout = ListenBeamform_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% --- Executes on selection change in lb_sourcesLeft.
function lb_sourcesLeft_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% If source changed, stop audio
stop(handles.audio1);


% Set list box to interruptible
% set(hObject, 'Interruptible', 'on');

% Get source(s) selected by user    source = [value]
sources = get(hObject,'Value');
handles.sources = sources;

% Update handles structure
guidata(hObject,handles);


% Hints: contents = cellstr(get(hObject,'String')) returns lb_sourcesLeft contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_sourcesLeft

%% --- Executes during object creation, after setting all properties.
function lb_sourcesLeft_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_sourcesLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in pb_playBF.
function pb_playBF_Callback(hObject, eventdata, handles)
% hObject    handle to pb_playBF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set default to Source1
if isempty(handles.sources)
    handles.sources = 1;
    disp('Stream set to default - Stream1')
end

% If audio is playing when 'PLAY' is selected, stop current file(s)
if isplaying(handles.audio1)
    stop(handles.audio1);
    disp('Current sound file(s) must be stopped before pressing play.');
end


% Audio data
y = handles.y;      fs = handles.fs;
[handles.rows cols] = size(y);

% Convert timestamps to samples to play with audioplayer
handles.startsamp(:,1) = floor(handles.StartTimes(:,1).*(fs(1,:)'));

numSources = length(handles.sources);       % Number of sources selected by user

% User selected source(s)
if numSources>1                        % User selected multiple sources
    % Find selected sources
    source = zeros(numSources,1);
    for n=1:numSources
        source(n,1) = handles.sources(1,n);
    end

    % Find smallest sample to start playing
    start = handles.startsamp(:,1);
    startnew = zeros(numSources,1);
    fsnew = zeros(numSources,1);
    for k=1:numSources
        value = source(k,1);
        startnew(k,1) = start(value,1);
        fsnew(k,1) = handles.fs(1,value);
    end
    
    % Check sampling rate for all sources chosen
    % If all are not the same, set to default of 44.1 kHz
    if (sum(fsnew)/numSources) ~= fsnew(1,1)            % If all are not the same
        handles.fsall(1,1) = 44100;                     % Use default FS of 44.1 kHz
    else
        handles.fsall(1,1) = fsnew(1,1);                % Else, all are same FS
    end

    % Find minimum value and set to start time
    beg = min(startnew);
    % Convert to samples
    beg = beg*handles.fsall;

    % Import audio for each selected source
    handles.x(:,1) = zeros(handles.maxY,1);
    for i=1:numSources
        value = source(i,1);
        handles.ynew(:,i) = handles.y(:,value);
        handles.fsnew(1,i) = handles.fsall(1,1);
        % Add streams together for chosen sources
        handles.x(:,1) = handles.x(:,1) + handles.ynew(:,i);
    end

    % Assign data to audioplayer object
    handles.audio1 = audioplayer(handles.x(:,1),fsnew(1,1));

    % If audio is playing when 'PLAY' is selected, stop current file(s)
    if isplaying(handles.audio1)
        stop(handles.audio1);
        disp('Current sound file(s) must be stopped before pressing play.');
    end
    
    % Play multiple audio
   % if beg == 0
   %     play(handles.audio1);
    %elseif beg ~= 0
        play(handles.audio1,beg);
   % else
    %    error('Error playing sound files!')
   % end
%*************************************************************************    
elseif numSources==1                    % User selected one source
    % Find selected sources
    source(1,1) = handles.sources(1,1);

    % Find smallest sample to start playing
    start = handles.startsamp(:,1);
    value = source(1,1);
    startnew(1,1) = start(value,1);
    fsnew(1,1) = handles.fs(1,value);
    handles.fsall(1,1) = fsnew(1,1);                
 
    % Find minimum value and set to start time
    beg = min(startnew);
    % Convert to samples
%     beg = beg*handles.fsall;    
    
    % Import audio for each selected source
    handles.ynew(:,1) = handles.y(:,value);
    handles.fsnew(1,1) = handles.fsall(1,1);
    
    % Assign sound data to each audioplayer object
    handles.audio1 = audioplayer(handles.ynew(:,1),fsnew(1,1));

    % If audio is playing when 'PLAY' is selected, stop current file(s)
    if isplaying(handles.audio1)
        stop(handles.audio1);
        disp('Current sound file(s) must be stopped before pressing play.');
    end
    
    
    % Play multiple audio
 %   if beg == 0
 %       play(handles.audio1);
 %   elseif beg ~= 0
        play(handles.audio1,beg);
  %  else
  %      error('Error playing sound file!')
  %  end
else
    error('Error playing sound file(s)!')
end


% Update handles
guidata(hObject,handles);


%% --- Executes on button press in pb_stopBF.
function pb_stopBF_Callback(hObject, eventdata, handles)
% hObject    handle to pb_stopBF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isplaying(handles.audio1)        % If playing, Stop file
    stop(handles.audio1);
elseif ~isplaying(handles.audio1)   % If not playing, display text
    disp('Sound file can only be stopped if previously playing.')
else
    disp('Error stopping sound file!')
end


guidata(hObject,handles);


% --- Executes on button press in pb_addBF.
function pb_addBF_Callback(hObject, eventdata, handles)
% hObject    handle to pb_addBF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
          
% Add stream to listbox on right, and remove from list box on left
allstreams = length(handles.leftstreams);       % Number of possible streams in left lb
leftStreams = handles.leftstreams;              % Streams in left lb
if isempty(handles.Maskstrings)                 % All streams in right lb
    handles.Maskstrings = [];
end

addstream = handles.sources;            % Stream selected to be added to rb
rmvstream = handles.sources;            % Stream selected to be removed from lb (same as above) - this is an index value to the array
maskTemp = leftStreams;
leftStreams(rmvstream) = 0;             % Set stream to be removed in lb = 0 using index
countR = 1;
countL = 1;
tempRight = [];
tempLeft = [];
for i=1:allstreams
    if i==addstream                     % If i equals stream selected to be added
        tempRight(1,countR)= maskTemp(i);          % Store i in temp var to add to right lb
        countR = countR+1;
    end
    if (leftStreams(i)~=0)              % Remove stream that is labeled 0
        tempLeft(1,countL) = leftStreams(i);
        countL = countL+1;
    end    
end


% Right List Box
handles.maskstreams = [handles.maskstreams, tempRight];     % Append tempR to Masked streams in right LB
handles.maskstreams = sort(handles.maskstreams);          % Sort in ascending order
handles.Maskstrings = [];               % Holds string names      
for st=1:length(handles.maskstreams)    % Convert to strings
    stream = handles.maskstreams(1,st);
    if stream<10
        handles.Maskstrings =  [ handles.Maskstrings; 'Stream0', num2str(stream)];
    elseif stream<100
        handles.Maskstrings =  [ handles.Maskstrings; 'Stream', num2str(stream)];
    end          
end
% Update right listbox to only have non deleted streams
set(handles.lb_maskstreams, 'String', handles.Maskstrings);     % Update right lb
set(handles.lb_maskstreams, 'Value' , 1);


% Left List Box
if isempty(tempLeft)                    % If no streams left in left LB
    handles.leftstreams = [];
else
    handles.leftstreams = tempLeft;     % Keep all streams saved in tempLeft
end
handles.leftstreams = sort(handles.leftstreams);          % Sort in ascending order
handles.Leftstrings = [];               % Holds string names 
for st=1:length(handles.leftstreams)    % Convert to strings
    stream = handles.leftstreams(1,st);
    if stream<10
        handles.Leftstrings =  [ handles.Leftstrings; 'Stream0', num2str(stream)];
    elseif stream<100
        handles.Leftstrings =  [ handles.Leftstrings; 'Stream', num2str(stream)];
    end          
end
% Update left listbox to contain the stream deleted from the right LB
set(handles.lb_sourcesLeft, 'String', handles.Leftstrings);       % Update left lb
set(handles.lb_sourcesLeft, 'Value' , 1);     % Set Value to 1 so listbox is not deleted

% Update Handles
guidata(hObject,handles);



%% --- Executes on button press in pb_masking.
function pb_masking_Callback(hObject, eventdata, handles)
% hObject    handle to pb_masking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf);
if(ishandle(handles.hSources)) 
close(handles.hSources);
end

% If no streams are selected to mask, display error and GUI again
if isempty(handles.maskstreams)
    disp('No streams were selected for masking! Select stream(s) to mask.')
    ListenBeamform({handles.bffname});          % Display GUI again
else
    moviemake_flag = handles.movieflag;
    trez = handles.trez;
    stseg = handles.stseg;
    fovp = handles.fovp;

    % Read in sources to do masking
    validsources = length(handles.maskstreams);

    % Delete beamformed wave files that are not listed
    for num=1:length(handles.allstreams)
        a = find(handles.allstreams(1,num) == handles.maskstreams(1,:));
        if isempty(a)        
            if num<10
                fntest = [handles.pn, handles.basenam, '0' int2str(num),'.wav'];     % Stream file name
            else
                fntest = [handles.pn, handles.basenam,  int2str(num),'.wav'];        % Stream file name
            end
            fileexist = exist(fntest,'file');
            if fileexist==2
                delete(fntest);
            end
        end
    end

    % Create masked streams and save in separate files
    scriptm = masksegment(handles.basenam, handles.pn, handles.scriptout, handles.maskstreams, validsources);

    %  Create graphical output ?
    if moviemake_flag == 1
        %  Create avi file to for movie of detected streams (must create a unique
        %  avi file name each time.)  Will not overwrite exisiting AVI file
        fratein = 1/(trez/2);                       % Frames per second on input
        frateout = 15;                              % Output frame rate
        fname = [handles.basenam '.avi'];           % File name of movie
        sfactor = 200;                              % Scale detection stat value to marker size
        mv = mkmov3d(stseg,fratein,frateout,fname,sfactor);
    end

    % Save scripting files for streams
    scfname = [handles.basenam 'scripts.mat'];
    save(scfname,'fovp','stseg','scriptm')          % Save file as .mat

    % Pass data to GUI to display and listen to sources
    ListenMult({scfname});
end


%% --------------------------------------------------------------------
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(gcf || handles.hSources)
    disp('Please close figures by clicking the ''X'' in the top right corner.')
else
    close(gcf);
    close(handles.hSources);
end


%% --- Executes on selection change in lb_maskstreams.
function lb_maskstreams_Callback(hObject, eventdata, handles)
% hObject    handle to lb_maskstreams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_maskstreams contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_maskstreams

% Stop audio first
% if isplaying(handles.audio1)
%     stop(handles.audio1);
% end

% Get source(s) selected by user
% If one source selected, sources = [value]
% If two or more sources selected, sources = [value1 value2 ...]
rmv = get(hObject,'Value');
handles.sourcermv = rmv;


% Update handles structure
guidata(hObject,handles);





%% --- Executes during object creation, after setting all properties.
function lb_maskstreams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_maskstreams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%% --- Executes on button press in pb_removeBF.
function pb_removeBF_Callback(hObject, eventdata, handles)
% hObject    handle to pb_removeBF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% If no kept streams to remove, error
if isempty(handles.maskstreams)             
    error('You must select a stream to remove from the listbox')
end

% Remove stream from right list box and add to left list box
totalmskstr = length(handles.maskstreams);  % Number of streams in right LB
totalmask = handles.maskstreams;            % Streams in left LB
if isempty(handles.Leftstrings)             % If no streams shown in left LB
    handles.Leftstrings = [];
end

rmvstream = handles.sourcermv;        % Stream selected to be removed from right LB - index value
addstream = handles.sourcermv;        % Stream selected to be added to left LB (same as above) - index value
leftTemp = totalmask;                 % Copy all streams in right LB
totalmask(rmvstream) = 0;             % Set stream to be removed in right LB to 0

countR = 1;
countL = 1;
tempRight =[];
tempLeft = [];
for i=1:totalmskstr
    if (totalmask(i)~=0)               % Remove stream from right LB
        tempRight(1,countR) = totalmask(i);
        countR =countR+1;
    end
    if i==addstream                        % If i equals stream selected to be added to left LB
        tempLeft(1,countL)= leftTemp(i);   % Store i in temp var to add to left lb
        countL = countL+1;
    end
end


% Right List Box
if isempty(tempRight)
    handles.maskstreams = [];
else
   handles.maskstreams = [tempRight];     % Keep all streams saved in tempRight
end
handles.maskstreams = sort(handles.maskstreams);    % Sort in ascending order
handles.Maskstrings = [];               % Holds string names      
for st=1:length(handles.maskstreams)    % Convert to strings
    stream = handles.maskstreams(1,st);
    if stream<10
        handles.Maskstrings =  [ handles.Maskstrings; 'Stream0', num2str(stream)];
    elseif stream<100
        handles.Maskstrings =  [ handles.Maskstrings; 'Stream', num2str(stream)];
    end          
end
% Update right listbox to only have non deleted streams
set(handles.lb_maskstreams, 'String', handles.Maskstrings);     % Update right lb
set(handles.lb_maskstreams, 'Value' , 1);
 

% Left List Box
handles.leftstreams = [handles.leftstreams tempLeft];     % Keep all streams saved in tempLeft
handles.leftstreams = sort(handles.leftstreams);          % Sort in ascending order
handles.Leftstrings = [];               % Holds string names       
for st=1:length(handles.leftstreams)    % Convert to strings
    stream = handles.leftstreams(1,st);
    if stream<10
        handles.Leftstrings =  [ handles.Leftstrings; 'Stream0', num2str(stream)];
    elseif stream<100
        handles.Leftstrings =  [ handles.Leftstrings; 'Stream', num2str(stream)];
    end          
end
% Update left listbox to contain the stream deleted from the right LB
set(handles.lb_sourcesLeft, 'String', handles.Leftstrings);       % Update left lb
set(handles.lb_sourcesLeft, 'Value' , 1);     % Set Value to 1 so listbox is not deleted


% Update Handles
guidata(hObject,handles);


