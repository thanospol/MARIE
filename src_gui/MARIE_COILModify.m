function varargout = MARIE_COILModify(varargin)
% MARIE_COILMODIFY MATLAB code for MARIE_COILModify.fig
%      MARIE_COILMODIFY, by itself, creates a new MARIE_COILMODIFY or raises the existing
%      singleton*.
%
%      H = MARIE_COILMODIFY returns the handle to a new MARIE_COILMODIFY or the handle to
%      the existing singleton*.
%
%      MARIE_COILMODIFY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_COILMODIFY.M with the given input arguments.
%
%      MARIE_COILMODIFY('Property','Value',...) creates a new MARIE_COILMODIFY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_COILModify_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_COILModify_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_COILModify

% Last Modified by GUIDE v2.5 11-Nov-2014 17:25:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_COILModify_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_COILModify_OutputFcn, ...
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


% --- Executes just before MARIE_COILModify is made visible.
function MARIE_COILModify_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_COILModify (see VARARGIN)

global COIL;
handles.axes = varargin{1}; % input parameter

% Choose default command line output for MARIE_COILModify
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_COILModify wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_COILModify_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function COIL_ModelName_Callback(hObject, eventdata, handles)
% hObject    handle to COIL_ModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of COIL_ModelName as text
%        str2double(get(hObject,'String')) returns contents of COIL_ModelName as a double


% --- Executes during object creation, after setting all properties.
function COIL_ModelName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to COIL_ModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global COIL;
if ~isempty(COIL.name)
    set(hObject,'String',COIL.name);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Apply.
function Apply_Callback(hObject, eventdata, handles)
% hObject    handle to Apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global COIL;

% get data
angle = str2num(get(handles.COIL_Angle, 'String'));

trans = [ str2num(get(handles.COIL_Xtrans, 'String')),...
    str2num(get(handles.COIL_Ytrans, 'String')),...
    str2num(get(handles.COIL_Ztrans, 'String'))];

[COIL] = Modify_COIL(COIL,angle,trans);
COIL.name = get(handles.COIL_ModelName, 'String');

      
uiwait(msgbox('Changes Applied','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function COIL_Xtrans_Callback(hObject, eventdata, handles)
% hObject    handle to COIL_Xtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of COIL_Xtrans as text
%        str2double(get(hObject,'String')) returns contents of COIL_Xtrans as a double


% --- Executes during object creation, after setting all properties.
function COIL_Xtrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to COIL_Xtrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0');


function COIL_Ytrans_Callback(hObject, eventdata, handles)
% hObject    handle to COIL_Ytrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of COIL_Ytrans as text
%        str2double(get(hObject,'String')) returns contents of COIL_Ytrans as a double


% --- Executes during object creation, after setting all properties.
function COIL_Ytrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to COIL_Ytrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0');

function COIL_Ztrans_Callback(hObject, eventdata, handles)
% hObject    handle to COIL_Ztrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of COIL_Ztrans as text
%        str2double(get(hObject,'String')) returns contents of COIL_Ztrans as a double


% --- Executes during object creation, after setting all properties.
function COIL_Ztrans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to COIL_Ztrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0');

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function COIL_Angle_Callback(hObject, eventdata, handles)
% hObject    handle to COIL_Angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of COIL_Angle as text
%        str2double(get(hObject,'String')) returns contents of COIL_Angle as a double


% --- Executes during object creation, after setting all properties.
function COIL_Angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to COIL_Angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String','0');

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
