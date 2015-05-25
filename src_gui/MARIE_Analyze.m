function varargout = MARIE_Analyze(varargin)
% MARIE_ANALYZE MATLAB code for MARIE_Analyze.fig
%      MARIE_ANALYZE, by itself, creates a new MARIE_ANALYZE or raises the existing
%      singleton*.
%
%      H = MARIE_ANALYZE returns the handle to a new MARIE_ANALYZE or the handle to
%      the existing singleton*.
%
%      MARIE_ANALYZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_ANALYZE.M with the given input arguments.
%
%      MARIE_ANALYZE('Property','Value',...) creates a new MARIE_ANALYZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_Analyze_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_Analyze_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_Analyze

% Last Modified by GUIDE v2.5 25-Jun-2014 19:46:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_Analyze_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_Analyze_OutputFcn, ...
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


% --- Executes just before MARIE_Analyze is made visible.
function MARIE_Analyze_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_Analyze (see VARARGIN)

% Choose default command line output for MARIE_Analyze
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_Analyze wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_Analyze_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
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
set(hObject,'String', RHBM.name);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


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
set(hObject,'String', COIL.name);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Close_Button.
function Close_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Close_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in FullMR_Button.
function FullMR_Button_Callback(hObject, eventdata, handles)
% hObject    handle to FullMR_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_SolveMR;
uiwait(gcf);

uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in COIL_FreqSweep_Button.
function COIL_FreqSweep_Button_Callback(hObject, eventdata, handles)
% hObject    handle to COIL_FreqSweep_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_SolveFreq;
uiwait(gcf);

uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in RHBM_Scat_Button.
function RHBM_Scat_Button_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Scat_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_SolveScat;
uiwait(gcf);

uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
