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

% Last Modified by GUIDE v2.5 29-Dec-2018 06:14:12

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
  handles.isStart = 0;
  handles.isRun = 0;
  handles.capture = 0;
  handles.isEnroll = 1;
  count = 0;
  time = 0;
  set(handles.buttonStop,'Enable','on');
  set(handles.buttonEnroll,'Enable','off');
  set(handles.buttonStart,'Enable','on');
  set(handles.captureButton,'Enable','on');
  set(handles.labelEdit,'Enable','on');
  set(handles.popupmenu,'Enable','on');
  guidata(hObject,handles);
  
  while (handles.isEnroll)
      tic();
      [~,img] = openCamera(handles);
      handles = guidata(hObject);  
      handles.img_real = img;
      guidata(hObject,handles);
      time = time + toc();
      count = count + 1;
      text(10, 10, ['FPS: ' int2str(floor(count/time))], 'color', [0 1 1]);
      bbox = py.facerecog.detect(img);
      bbox = cell(bbox);
      if py.len(bbox) > 0
        if py.len(bbox) == 1
            set(handles.captureButton,'Enable','on');
        else
            set(handles.captureButton,'Enable','off');
        end
        for n = 1:py.len(bbox)
            top = double(bbox{n}{1})*4;
            right = double(bbox{n}{2})*4;
            bottom = double(bbox{n}{3})*4;
            left = double(bbox{n}{4})*4;
            b = [left, top, right-left ,bottom-top];
            rect = rectangle('Position',b);
            set(rect,'FaceColor','none','EdgeColor','b','LineWidth',2);
            rect_title(rect,"Face");
            
            if (handles.capture)
                folder = 'images/';
                label = handles.labelEdit.String;
                fileType = '.jpg';
                s = strcat(folder,label,fileType);
                img = handles.img_real;
                img = imcrop(img,b);
                imwrite(img,s);
                imshow(img)
                handles.capture = 0;
                handles.isEnroll = 0;
                break;
            end
                       
        end
      else
          set(handles.captureButton,'Enable','off');
      end
      if (time >= 1)
        time = 0;
        count = 0;
      end
      axis image;
  end
  
  guidata(hObject, handles); % Update handles structure

% --- Executes on button press in buttonStart.
function buttonStart_Callback(hObject, eventdata, handles)
% hObject    handle to buttonStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setup_paths();
handles.isStart = 1;
handles.isEnroll = 0;
count = 0;
time = 0;
handles.isRun = 0;
handles.isFound = 0;
handles.countList = 0;
handles.listParams = struct([]);
handles.drawBbox = [];
set(handles.buttonStop,'Enable','on');
set(handles.buttonEnroll,'Enable','on');
set(handles.buttonStart,'Enable','off');
set(handles.captureButton,'Enable','off');
set(handles.labelEdit,'Enable','off');
set(handles.popupmenu,'Enable','off');
%system('python faceregcog.py');
py.facerecog.init();
guidata(hObject,handles);

while (handles.isStart)
    tic();
    %drawnow %Give the button callback a chance to interrupt the 
    [img_gray,img_real] = openCamera(handles);
    handles = guidata(hObject);  %Get the newest GUI data
    handles.img_real = img_real;
    bbox = py.facerecog.recog(img_real);
    bbox = cell(bbox);
    if py.len(bbox) > 0
        for n = 1:py.len(bbox)
            name = string(bbox{n}{2});
            if name ~= "Unknown" 
                top = double(bbox{n}{1}{1})*4;
                right = double(bbox{n}{1}{2})*4;
                bottom = double(bbox{n}{1}{3})*4;
                left = double(bbox{n}{1}{4})*4;
                b = [left, top, right-left ,bottom-top];              
                indexSame = -1;
                [~,handles.countList] = size(handles.listParams);
                for i = 1:handles.countList
                    if handles.listParams(i).params.label == name
                        indexSame = i;
                        break;
                    end
                end
                if indexSame ~= -1
                    handles.listParams(indexSame) = [];;
                end
                [params] = create_realtime_sequence(img_gray,b);
                params.color = "g";
                params.label = name;
                handles.params = run_realtime_STRCF(params, handles);
                [~,handles.countList] = size(handles.listParams);
                handles.listParams(handles.countList+1).params = handles.params;   
            end
        end
    end
    [~,handles.countList] = size(handles.listParams);
    if handles.countList > 0 && ~isempty(handles.countList)
        indexOut = -1;       
        for i = 1:handles.countList
            handles.listParams(i).params.im = img_gray;
            handles.listParams(i).params = run_realtime_STRCF(handles.listParams(i).params,handles);           
            bbox = handles.listParams(i).params.rect_position_vis;
            [height, ~] = size(handles.img_real);
            if (0 - bbox(1,3)/2 > bbox(1,1) || 0 - bbox(1,4)/2 > bbox(1,2) || bbox(1,1) + bbox(1,3)/2 > 640 || bbox(1,2) + bbox(1,4)/2 > height )
                indexOut = i;
            end
            rect = rectangle('Position',handles.listParams(i).params.rect_position_vis);
            set(rect,'FaceColor','none','EdgeColor',handles.listParams(i).params.color,'LineWidth',2);
            rect_title(rect,handles.listParams(i).params.label);
        end
        if indexOut ~= -1
            handles.listParams(indexOut) = [];            
            [~,handles.countList] = size(handles.listParams);            
            guidata(hObject,handles);
        end
    end
    time = time + toc();
    count = count + 1;
    text(10, 10, ['FPS: ' int2str(floor(count/time))], 'color', [0 1 1]);
    if (time >= 2)
        time = 0;
        count = 0;
    end
    axis image;
    guidata(hObject,handles)
   
end

function rect_title(rect,str)
global updateflag
updateflag = false;
% Add initial title to UserData
rect.UserData.str = str;
rect.UserData.t = text(rect.Position(1)+4,rect.Position(2)+21,rect.UserData.str,'HorizontalAlignment','left','FontSize',25, 'BackgroundColor', 'g');
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

function [img_gray,img] = openCamera(handles)
    axes(handles.image)
    
    img = snapshot(handles.cam);
    img = imresize(img,[480 640]);
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
  handles.isEnroll = 0;
  set(handles.buttonStop,'Enable','off');
  set(handles.buttonEnroll,'Enable','off');
  set(handles.captureButton,'Enable','off');
  set(handles.buttonStart,'Enable','on');
  set(handles.popupmenu,'Enable','off');
  set(handles.labelEdit,'Enable','off');
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


% --- Executes on button press in captureButton.
function captureButton_Callback(hObject, eventdata, handles)
% hObject    handle to captureButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.capture = 1;
guidata(hObject, handles);
