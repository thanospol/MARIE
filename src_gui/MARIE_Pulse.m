function varargout = MARIE_Pulse(varargin)
% MARIE_PULSE MATLAB code for MARIE_Pulse.fig
%      MARIE_PULSE, by itself, creates a new MARIE_PULSE or raises the existing
%      singleton*.
%
%      H = MARIE_PULSE returns the handle to a new MARIE_PULSE or the handle to
%      the existing singleton*.
%
%      MARIE_PULSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_PULSE.M with the given input arguments.
%
%      MARIE_PULSE('Property','Value',...) creates a new MARIE_PULSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_Pulse_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_Pulse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_Pulse

% Last Modified by GUIDE v2.5 22-Dec-2014 16:44:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_Pulse_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_Pulse_OutputFcn, ...
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


% --- Executes just before MARIE_Pulse is made visible.
function MARIE_Pulse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_Pulse (see VARARGIN)

% Choose default command line output for MARIE_Pulse
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_Pulse wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_Pulse_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Body_Model_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Body_Model_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Body_Model_Edit as text
%        str2double(get(hObject,'String')) returns contents of Body_Model_Edit as a double


% --- Executes during object creation, after setting all properties.
function Body_Model_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Body_Model_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if strcmp(RHBM.name, 'No Model Selected')
    set(hObject,'String','Free Space');
else
    set(hObject,'String',RHBM.name);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Coil_Model_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Coil_Model_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Coil_Model_Edit as text
%        str2double(get(hObject,'String')) returns contents of Coil_Model_Edit as a double


% --- Executes during object creation, after setting all properties.
function Coil_Model_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Coil_Model_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global COIL;
set(hObject,'String',COIL.name);

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
function Pulse_File_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Pulse_File_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pulse_File_Edit as text
%        str2double(get(hObject,'String')) returns contents of Pulse_File_Edit as a double


% --- Executes during object creation, after setting all properties.
function Pulse_File_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pulse_File_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Visualize_Button.
function Visualize_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Visualize_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_Visualize_pulse;
uiwait(gcf);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Apply_Button.
function Apply_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Apply_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RHBM;
global COIL;
global SOL;
global PULSE;

% get file name with the pulse sequence
filename = get(handles.Pulse_File_Edit, 'String');

% call the function to apply the pulse sequence
[solution_pulse] = pulse_apply(SOL,COIL,RHBM,filename);

if isempty(solution_pulse)
    f = warndlg('Pulse not applied', 'Pulse not applied');
    waitfor(f);
else
    PULSE = solution_pulse;
end

uiwait(msgbox('Pulse Sequence Applied','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
uiresume(gcbf);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Select_File_Pulse.
function Select_File_Pulse_Callback(hObject, eventdata, handles)
% hObject    handle to Select_File_Pulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get file name
[pulse_filename, pulse_path] = uigetfile('*.mps', 'Select file with Pulse Sequence');
% check that return was correct
if (pulse_filename(1) ~= 0) && (pulse_path(1) ~= 0)
    set(handles.Pulse_File_Edit, 'String', fullfile(pulse_path, pulse_filename));
else
    f = warndlg('Incorrect file selected', 'Wrong File');
    waitfor(f);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
