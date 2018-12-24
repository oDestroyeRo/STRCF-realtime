function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 22-Dec-2018 10:00:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)
handles.cam = webcam();

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function labelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to labelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of labelEdit as text
%        str2double(get(hObject,'String')) returns contents of labelEdit as a double


% --- Executes during object creation, after setting all properties.
function labelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in buttonEnroll.
function buttonEnroll_Callback(hObject, eventdata, handles)
% hObject    handle to buttonEnroll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch handles.popupmenu.Value
    case 1
        color = 'r';
    case 2
        color = 'g';
    case 3
        color = 'b';
    case 4
        color = 'y';               
end  
handles.drawBbox = drawrectangle(handles.image,'LineWidth',5,'Color',color);
handles.drawBbox = floor(handles.drawBbox.Position);
guidata(hObject,handles);

% --- Executes on button press in buttonStart.
function buttonStart_Callback(hObject, eventdata, handles)
% hObject    handle to buttonStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setup_paths();
handles.isStart = 1;
count = 1;
time = 0.1;
handles.isRun = 0;
handles.isFound = 0;
handles.countList = 0;
handles.listParams = struct([]);
handles.drawBbox = [];
set(handles.buttonStop,'Enable','on');
set(handles.buttonEnroll,'Enable','on');
set(handles.buttonStart,'Enable','off');
guidata(hObject,handles);

while (handles.isStart)
    tic();
    %drawnow %Give the button callback a chance to interrupt the 
    img_gray = openCamera(handles);
    handles = guidata(hObject);  %Get the newest GUI data
    
    if ~isempty(handles.drawBbox)
        [params] = create_realtime_sequence(img_gray,handles.drawBbox);
        switch handles.popupmenu.Value
            case 1
                params.color = 'r';
            case 2
                params.color = 'g';
            case 3
                params.color = 'b';
            case 4
                params.color = 'y';               
        end       
        params.label = handles.labelEdit.String;
        handles.params = run_realtime_STRCF(params, handles);
        [~,sizeParams] = size(handles.listParams);
        handles.countList = sizeParams;
        handles.listParams(handles.countList+1).params = handles.params;
        handles.drawBbox = [];
        handles.countList = sizeParams;
    end
    if ~isempty(handles.listParams)
        listParams = handles.listParams;
        countList = handles.countList;
        for i = 1:countList+1
            listParams(i).params.im = img_gray;
            listParams(i).params = run_realtime_STRCF(listParams(i).params,handles);
            rect = rectangle('Position',listParams(i).params.rect_position_vis);
            set(rect,'FaceColor','none','EdgeColor',listParams(i).params.color,'LineWidth',1);
            rect_title(rect,listParams(i).params.label);
        end
        handles.listParams = listParams;
        %rectangle('Position',handles.listParams(1).params.rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
    end 
    %listParams = handles.listParams;
    text(10, 10, ['FPS: ' int2str(count/time)], 'color', [0 1 1]);
    axis image;
    count = count + 1;
    time = time + toc();
    guidata(hObject,handles)
   
end

function rect_title(rect,str)
global updateflag
updateflag = false;
% Add initial title to UserData
rect.UserData.str = str;
rect.UserData.t = text(rect.Position(1)+10,rect.Position(2)-10,rect.UserData.str,'HorizontalAlignment','center');
% Add Callbacks
ax = rect.Parent;
f = ax.Parent;
f.WindowButtonUpFcn = @draw_rectangle;
rect.ButtonDownFcn = @draw_rectangle;

function draw_rectangle(src,evt)
global updateflag
if ishandle(src) && strcmp(get(src,'type'),'figure') && updateflag
    % Assign temporary data to rect
    RectData = src.UserData.RectData;
    % Move Rectangle
    Pold = RectData.rect.Position;
    P1 = RectData.Pdown.IntersectionPoint;
    P2 = evt.IntersectionPoint;
    dx = P2(1) - P1(1);
    dy = P2(2) - P1(2);
    RectData.rect.Position = Pold + [dx dy 0 0];
    % Update title
    delete(RectData.t)
    RectData.rect.UserData.t = text(RectData.rect.Position(1)+RectData.rect.Position(3)/2,RectData.rect.Position(2)+RectData.rect.Position(4)/2,RectData.rect.UserData.str,'HorizontalAlignment','center');
    updateflag = false;
elseif ishandle(src) && strcmp(get(src,'type'),'rectangle')
    updateflag = true;
    ax = src.Parent;
    f = ax.Parent;
    f.UserData.RectData.t = src.UserData.t;
    f.UserData.RectData.rect = src;
    f.UserData.RectData.Pdown = evt;
end

function img_gray = openCamera(handles)
    axes(handles.image)
    img = snapshot(handles.cam);
    img = flip(img,2);
    img_gray = rgb2gray(img);
    imagesc(img);
    
    %axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);   

% --- Executes on button press in buttonStop.
function buttonStop_Callback(hObject, eventdata, handles)
% hObject    handle to buttonStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  handles.isStart = 0;
  handles.isRun = 0;
  set(handles.buttonStop,'Enable','off');
  set(handles.buttonEnroll,'Enable','off');
  set(handles.buttonStart,'Enable','on');
  guidata(hObject, handles); % Update handles structure


% --- Executes on selection change in popupmenu.
function popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu


% --- Executes during object creation, after setting all properties.
function popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
