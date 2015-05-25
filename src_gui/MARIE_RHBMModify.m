function varargout = MARIE_RHBMModify(varargin)
% MARIE_RHBMMODIFY MATLAB code for MARIE_RHBMModify.fig
%      MARIE_RHBMMODIFY, by itself, creates a new MARIE_RHBMMODIFY or raises the existing
%      singleton*.
%
%      H = MARIE_RHBMMODIFY returns the handle to a new MARIE_RHBMMODIFY or the handle to
%      the existing singleton*.
%
%      MARIE_RHBMMODIFY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_RHBMMODIFY.M with the given input arguments.
%
%      MARIE_RHBMMODIFY('Property','Value',...) creates a new MARIE_RHBMMODIFY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_RHBMModify_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_RHBMModify_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_RHBMModify

% Last Modified by GUIDE v2.5 25-Jun-2014 19:53:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_RHBMModify_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_RHBMModify_OutputFcn, ...
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


% --- Executes just before MARIE_RHBMModify is made visible.
function MARIE_RHBMModify_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_RHBMModify (see VARARGIN)

global RHBM;
handles.axes = varargin{1}; % input parameter

% Choose default command line output for MARIE_RHBMModify
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_RHBMModify wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_RHBMModify_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function RHBM_ModelName_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_ModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_ModelName as text
%        str2double(get(hObject,'String')) returns contents of RHBM_ModelName as a double


% --- Executes during object creation, after setting all properties.
function RHBM_ModelName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_ModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.name)
    set(hObject,'String',RHBM.name);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Apply.
function Apply_Callback(hObject, eventdata, handles)
% hObject    handle to Apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RHBM;

% get data
resolution = [ str2num(get(handles.RHBM_Resolution, 'String'))];
xlimits = [ str2num(get(handles.RHBM_Xmin, 'String')),...
    str2num(get(handles.RHBM_Xmax, 'String'))];
ylimits = [ str2num(get(handles.RHBM_Ymin, 'String')),...
    str2num(get(handles.RHBM_Ymax, 'String'))];
zlimits = [ str2num(get(handles.RHBM_Zmin, 'String')),...
    str2num(get(handles.RHBM_Zmax, 'String'))];
trans = [ str2num(get(handles.RHBM_Xtrans, 'String')),...
    str2num(get(handles.RHBM_Ytrans, 'String')),...
    str2num(get(handles.RHBM_Ztrans, 'String'))];

[RHBM] = Resize_RHBM(RHBM,resolution,xlimits,ylimits,zlimits,trans);
RHBM.name = get(handles.RHBM_ModelName, 'String');

      
uiwait(msgbox('Changes Applied','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


function RHBM_Xtrans_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Xtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Xtrans as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Xtrans as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Xtrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Xtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0');


function RHBM_Ytrans_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Ytrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Ytrans as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Ytrans as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Ytrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Ytrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0');

function RHBM_Ztrans_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Ztrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Ztrans as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Ztrans as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Ztrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Ztrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0');


function RHBM_Resolution_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Resolution as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Resolution as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.r)
    x = RHBM.r(:,:,:,1);
    res = x(2,1,1) - x(1,1,1);
    value = sprintf('%.4g', res);
    set(hObject,'String',value);
end



function RHBM_Xmin_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Xmin as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Xmin as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Xmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Xmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.r)
    x = RHBM.r(:,:,:,1);
    xmin = min(x(:));
    value = sprintf('%.3g', xmin);
    set(hObject,'String',value);
end


function RHBM_Ymin_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Ymin as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Ymin as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Ymin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Ymin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.r)
    y = RHBM.r(:,:,:,2);
    ymin = min(y(:));
    value = sprintf('%.3g', ymin);
    set(hObject,'String',value);
end


function RHBM_Zmin_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Zmin as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Zmin as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Zmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Zmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.r)
    z = RHBM.r(:,:,:,3);
    zmin = min(z(:));
    value = sprintf('%.3g', zmin);
    set(hObject,'String',value);
end


function RHBM_Xmax_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Xmax as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Xmax as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Xmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Xmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.r)
    x = RHBM.r(:,:,:,1);
    xmax = max(x(:));
    value = sprintf('%.3g', xmax);
    set(hObject,'String',value);
end


function RHBM_Ymax_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Ymax as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Ymax as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Ymax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Ymax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.r)
    y = RHBM.r(:,:,:,2);
    ymax = max(y(:));
    value = sprintf('%.3g', ymax);
    set(hObject,'String',value);
end


function RHBM_Zmax_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Zmax as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Zmax as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Zmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Zmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.r)
    z = RHBM.r(:,:,:,3);
    zmax = max(z(:));
    value = sprintf('%.3g', zmax);
    set(hObject,'String',value);
end


% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcbf);
close(handles.figure1);
