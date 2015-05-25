function varargout = MARIE_SolveScat(varargin)
% MARIE_SOLVESCAT MATLAB code for MARIE_SolveScat.fig
%      MARIE_SOLVESCAT, by itself, creates a new MARIE_SOLVESCAT or raises the existing
%      singleton*.
%
%      H = MARIE_SOLVESCAT returns the handle to a new MARIE_SOLVESCAT or the handle to
%      the existing singleton*.
%
%      MARIE_SOLVESCAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_SOLVESCAT.M with the given input arguments.
%
%      MARIE_SOLVESCAT('Property','Value',...) creates a new MARIE_SOLVESCAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_SolveScat_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_SolveScat_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_SolveScat

% Last Modified by GUIDE v2.5 25-Jun-2014 21:27:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_SolveScat_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_SolveScat_OutputFcn, ...
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


% --- Executes just before MARIE_SolveScat is made visible.
function MARIE_SolveScat_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_SolveScat (see VARARGIN)

% Choose default command line output for MARIE_SolveScat
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_SolveScat wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_SolveScat_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function SCAT_ModelName_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_ModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_ModelName as text
%        str2double(get(hObject,'String')) returns contents of SCAT_ModelName as a double


% --- Executes during object creation, after setting all properties.
function SCAT_ModelName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_ModelName (see GCBO)
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
function Freq_Min_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_Min_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Freq_Min_Edit as text
%        str2double(get(hObject,'String')) returns contents of Freq_Min_Edit as a double


% --- Executes during object creation, after setting all properties.
function Freq_Min_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq_Min_Edit (see GCBO)
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
function Freq_Max_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_Max_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Freq_Max_Edit as text
%        str2double(get(hObject,'String')) returns contents of Freq_Max_Edit as a double


% --- Executes during object creation, after setting all properties.
function Freq_Max_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq_Max_Edit (see GCBO)
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
function Freq_Num_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_Num_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Freq_Num_Edit as text
%        str2double(get(hObject,'String')) returns contents of Freq_Num_Edit as a double


% --- Executes during object creation, after setting all properties.
function Freq_Num_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq_Num_Edit (see GCBO)
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
% --- Executes on button press in Freq_Linscale_Select.
function Freq_Linscale_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_Linscale_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Freq_Linscale_Select

set(handles.Freq_Logscale_Select, 'Value', 0);

% --- Executes on button press in Freq_Logscale_Select.
function Freq_Logscale_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_Logscale_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Freq_Logscale_Select

set(handles.Freq_Linscale_Select, 'Value', 0);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Tolerance_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Tolerance_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tolerance_Box as text
%        str2double(get(hObject,'String')) returns contents of Tolerance_Box as a double


% --- Executes during object creation, after setting all properties.
function Tolerance_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tolerance_Box (see GCBO)
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
% --- Executes on button press in Efield_check.
function Efield_check_Callback(hObject, eventdata, handles)
% hObject    handle to Efield_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Efield_check


% --- Executes on button press in Bfield_check.
function Bfield_check_Callback(hObject, eventdata, handles)
% hObject    handle to Bfield_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Bfield_check


% --- Executes on button press in Power_check.
function Power_check_Callback(hObject, eventdata, handles)
% hObject    handle to Power_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Power_check

% --- Executes on button press in SAR_check.
function SAR_check_Callback(hObject, eventdata, handles)
% hObject    handle to SAR_check (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SAR_check

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Dipole_Select.
function Dipole_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Dipole_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Dipole_Select

set(handles.PW_Select, 'Value', 0);
set(handles.CCL_Select, 'Value', 0);

set(handles.SCAT_PlaneWaveZpol, 'String', ' ');
set(handles.SCAT_PlaneWaveYpol, 'String', ' ');
set(handles.SCAT_PlaneWaveXpol, 'String', ' ');
set(handles.SCAT_PlaneWaveZdir, 'String', ' ');
set(handles.SCAT_PlaneWaveYdir, 'String', ' ');
set(handles.SCAT_PlaneWaveXdir, 'String', ' ');

set(handles.SCAT_DipoleZpos, 'String', '1.0');
set(handles.SCAT_DipoleYpos, 'String', '1.0');
set(handles.SCAT_DipoleXpos, 'String', '1.0');
set(handles.SCAT_DipoleZcomp, 'String', '1.0');
set(handles.SCAT_DipoleYcomp, 'String', '0.0');
set(handles.SCAT_DipoleXcomp, 'String', '0.0');

set(handles.CCL_Filename_Box, 'String', ' ');
set(handles.CCL_Amplitude_Box, 'String', ' ');

% --- Executes on button press in PW_Select.
function PW_Select_Callback(hObject, eventdata, handles)
% hObject    handle to PW_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PW_Select

set(handles.Dipole_Select, 'Value', 0);
set(handles.CCL_Select, 'Value', 0);

set(handles.SCAT_PlaneWaveZpol, 'String', '0.0');
set(handles.SCAT_PlaneWaveYpol, 'String', '0.0');
set(handles.SCAT_PlaneWaveXpol, 'String', '1.0');
set(handles.SCAT_PlaneWaveZdir, 'String', '1.0');
set(handles.SCAT_PlaneWaveYdir, 'String', '0.0');
set(handles.SCAT_PlaneWaveXdir, 'String', '0.0');

set(handles.SCAT_DipoleZpos, 'String', ' ');
set(handles.SCAT_DipoleYpos, 'String', ' ');
set(handles.SCAT_DipoleXpos, 'String', ' ');
set(handles.SCAT_DipoleZcomp, 'String', ' ');
set(handles.SCAT_DipoleYcomp, 'String', ' ');
set(handles.SCAT_DipoleXcomp, 'String', ' ');

set(handles.CCL_Filename_Box, 'String', ' ');
set(handles.CCL_Amplitude_Box, 'String', ' ');


% --- Executes on button press in CCL_Select.
function CCL_Select_Callback(hObject, eventdata, handles)
% hObject    handle to CCL_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CCL_Select

set(handles.PW_Select, 'Value', 0);
set(handles.Dipole_Select, 'Value', 0);

set(handles.SCAT_PlaneWaveZpol, 'String', ' ');
set(handles.SCAT_PlaneWaveYpol, 'String', ' ');
set(handles.SCAT_PlaneWaveXpol, 'String', ' ');
set(handles.SCAT_PlaneWaveZdir, 'String', ' ');
set(handles.SCAT_PlaneWaveYdir, 'String', ' ');
set(handles.SCAT_PlaneWaveXdir, 'String', ' ');

set(handles.SCAT_DipoleZpos, 'String', ' ');
set(handles.SCAT_DipoleYpos, 'String', ' ');
set(handles.SCAT_DipoleXpos, 'String', ' ');
set(handles.SCAT_DipoleZcomp, 'String', ' ');
set(handles.SCAT_DipoleYcomp, 'String', ' ');
set(handles.SCAT_DipoleXcomp, 'String', ' ');

% get file name
[ccl_filename, ccl_path] = uigetfile('*.wsd', 'Select file with Loop model (wire format)');
% check that return was correct
if (ccl_filename(1) ~= 0) && (ccl_path(1) ~= 0)
    set(handles.CCL_Filename_Box, 'String', fullfile(ccl_path, ccl_filename));
else
    f = warndlg('Incorrect file selected', 'Wrong File');
    waitfor(f);
end
set(handles.CCL_Amplitude_Box, 'String', '1.0');

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function SCAT_PlaneWaveZpol_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveZpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_PlaneWaveZpol as text
%        str2double(get(hObject,'String')) returns contents of SCAT_PlaneWaveZpol as a double


% --- Executes during object creation, after setting all properties.
function SCAT_PlaneWaveZpol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveZpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SCAT_PlaneWaveZdir_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveZdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_PlaneWaveZdir as text
%        str2double(get(hObject,'String')) returns contents of SCAT_PlaneWaveZdir as a double


% --- Executes during object creation, after setting all properties.
function SCAT_PlaneWaveZdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveZdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SCAT_PlaneWaveYpol_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveYpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_PlaneWaveYpol as text
%        str2double(get(hObject,'String')) returns contents of SCAT_PlaneWaveYpol as a double


% --- Executes during object creation, after setting all properties.
function SCAT_PlaneWaveYpol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveYpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SCAT_PlaneWaveYdir_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveYdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_PlaneWaveYdir as text
%        str2double(get(hObject,'String')) returns contents of SCAT_PlaneWaveYdir as a double


% --- Executes during object creation, after setting all properties.
function SCAT_PlaneWaveYdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveYdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SCAT_PlaneWaveXpol_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveXpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_PlaneWaveXpol as text
%        str2double(get(hObject,'String')) returns contents of SCAT_PlaneWaveXpol as a double


% --- Executes during object creation, after setting all properties.
function SCAT_PlaneWaveXpol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveXpol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SCAT_PlaneWaveXdir_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveXdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_PlaneWaveXdir as text
%        str2double(get(hObject,'String')) returns contents of SCAT_PlaneWaveXdir as a double


% --- Executes during object creation, after setting all properties.
function SCAT_PlaneWaveXdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_PlaneWaveXdir (see GCBO)
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
function SCAT_DipoleXpos_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleXpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_DipoleXpos as text
%        str2double(get(hObject,'String')) returns contents of SCAT_DipoleXpos as a double


% --- Executes during object creation, after setting all properties.
function SCAT_DipoleXpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleXpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SCAT_DipoleXcomp_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleXcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_DipoleXcomp as text
%        str2double(get(hObject,'String')) returns contents of SCAT_DipoleXcomp as a double


% --- Executes during object creation, after setting all properties.
function SCAT_DipoleXcomp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleXcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function SCAT_DipoleYpos_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleYpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_DipoleYpos as text
%        str2double(get(hObject,'String')) returns contents of SCAT_DipoleYpos as a double


% --- Executes during object creation, after setting all properties.
function SCAT_DipoleYpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleYpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SCAT_DipoleYcomp_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleYcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_DipoleYcomp as text
%        str2double(get(hObject,'String')) returns contents of SCAT_DipoleYcomp as a double


% --- Executes during object creation, after setting all properties.
function SCAT_DipoleYcomp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleYcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SCAT_DipoleZpos_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleZpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_DipoleZpos as text
%        str2double(get(hObject,'String')) returns contents of SCAT_DipoleZpos as a double


% --- Executes during object creation, after setting all properties.
function SCAT_DipoleZpos_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleZpos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SCAT_DipoleZcomp_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleZcomp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SCAT_DipoleZcomp as text
%        str2double(get(hObject,'String')) returns contents of SCAT_DipoleZcomp as a double


% --- Executes during object creation, after setting all properties.
function SCAT_DipoleZcomp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SCAT_DipoleZcomp (see GCBO)
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
function CCL_Filename_Box_Callback(hObject, eventdata, handles)
% hObject    handle to CCL_Filename_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CCL_Filename_Box as text
%        str2double(get(hObject,'String')) returns contents of CCL_Filename_Box as a double


% --- Executes during object creation, after setting all properties.
function CCL_Filename_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CCL_Filename_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function CCL_Amplitude_Box_Callback(hObject, eventdata, handles)
% hObject    handle to CCL_Amplitude_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CCL_Amplitude_Box as text
%        str2double(get(hObject,'String')) returns contents of CCL_Amplitude_Box as a double


% --- Executes during object creation, after setting all properties.
function CCL_Amplitude_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CCL_Amplitude_Box (see GCBO)
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
% --- Executes on button press in SCAT_Solve_Button.
function SCAT_Solve_Button_Callback(hObject, eventdata, handles)
% hObject    handle to SCAT_Solve_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global RHBM;
global SOL;

[L,M,N,~] = size(RHBM.r);

% get frequency data
numfreqs = round(str2num(get(handles.Freq_Num_Edit, 'String')));
freqmin = str2num(get(handles.Freq_Min_Edit, 'String'));
if (isempty(freqmin))
    errordlg('Wrong Frequency Range -- Please correct', 'Error');
end
if (numfreqs > 1) % frequency sweep
    freqmax = str2num(get(handles.Freq_Max_Edit, 'String'));
    if (isempty(freqmax)) || (freqmax <= freqmin)
        errordlg('Wrong Frequency Range -- Please correct', 'Error');
    end
    
    if get(handles.Freq_Logscale_Select,'Value') % logscale
        if (freqmin == 0)
            freqmin = min( [freqmax/100, 1e-6]); % avoid log of zero
        end
        fvalues = logspace(log10(freqmin),log10(freqmax),numfreqs);

    else % linear scale
        fvalues = linspace(freqmin,freqmax,numfreqs);
    end
else
    fvalues = freqmin; % single frequency
end
        
% store the frequencies
SOL.freq = squeeze(fvalues);

% get tolerance
tol = str2num(get(handles.Tolerance_Box, 'String'));


% get kind of excitation
type = find([get(handles.PW_Select, 'Value') get(handles.Dipole_Select, 'Value') get(handles.CCL_Select, 'Value')]);

switch type

case 1
    Exctype = 'P'; % plane wave excitation

    % direction
    vec1 = [ str2num(get(handles.SCAT_PlaneWaveXdir, 'String')),...
        str2num(get(handles.SCAT_PlaneWaveYdir, 'String')),...
        str2num(get(handles.SCAT_PlaneWaveZdir, 'String'))];

    % polarization
    vec2 = [ str2num(get(handles.SCAT_PlaneWaveXpol, 'String')),...
        str2num(get(handles.SCAT_PlaneWaveYpol, 'String')),...
        str2num(get(handles.SCAT_PlaneWaveZpol, 'String'))];
    

case 2
    Exctype = 'D'; % dipole excitation

    % position
    vec1 = [ str2num(get(handles.SCAT_DipoleXpos, 'String')),...
        str2num(get(handles.SCAT_DipoleYpos, 'String')),...
        str2num(get(handles.SCAT_DipoleZpos, 'String'))];

    % components
    vec2 = [ str2num(get(handles.SCAT_DipoleXcomp, 'String')),...
        str2num(get(handles.SCAT_DipoleYcomp, 'String')),...
        str2num(get(handles.SCAT_DipoleZcomp, 'String'))];

case 3
    Exctype = 'L'; % loop excitation

    % file
    vec1 = get(handles.CCL_Filename_Box, 'String');
    % amplitude
    vec2 = str2num(get(handles.CCL_Amplitude_Box, 'String'));


otherwise
    Exctype = [];

end


% generate excitation for multiple frequencies and solve
if ~isempty(Exctype)

    SOL.Jsol = zeros(L,M,N,3,1,numfreqs);
    
    if get(handles.Efield_check, 'Value')
        SOL.Esol = zeros(L,M,N,3,1,numfreqs);
    end
    if get(handles.Bfield_check, 'Value')
        SOL.Bsol = zeros(L,M,N,3,1,numfreqs);
    end
    if get(handles.SAR_check, 'Value')
        SOL.Ssol = zeros(L,M,N,1,numfreqs);
        SOL.Gsar = zeros(1,numfreqs);
    end
    if get(handles.Power_check, 'Value')
        SOL.Pabs = zeros(1,numfreqs);
    end

    for ii = 1:numfreqs
        % -------------------------------------------------
        % excitation
        [Einc,Hinc] = Generate_Excitation(RHBM.r,fvalues(ii),Exctype,vec1,vec2);
        fprintf(1,'\n\n ---------------------------------------------------------------------');
        fprintf(1,'\n ---------------------------------------------------------------------');
        fprintf(1,'\n Solving the VIE system at freq %.3f MHz\n ',fvalues(ii)/1e6);        
        % -------------------------------------------------
        %   Generate circulants if needed
        if (isempty(RHBM.fN))
            % we need to compute a new circulant
            [RHBM.fN] = getOPERATORS(RHBM.r,fvalues(ii),'N');
            RHBM.freqfN = fvalues(ii);
        else
            if (RHBM.freqfN ~= fvalues(ii))
                [RHBM.fN] = getOPERATORS(RHBM.r,fvalues(ii),'N');
                RHBM.freqfN = fvalues(ii);
            end
        end
        if (isempty(RHBM.fK))
            % we need to compute a new circulant
            [RHBM.fK] = getOPERATORS(RHBM.r,fvalues(ii),'K');
            RHBM.freqfK = fvalues(ii);
        else
            if (RHBM.freqfK ~= fvalues(ii))
                [RHBM.fK] = getOPERATORS(RHBM.r,fvalues(ii),'K');
                RHBM.freqfK = fvalues(ii);
            end
        end        
        % Solve
        [Jout,Sout,Eout,Bout,Gsar,Pabs,~,~] = Scat_Solver(Einc,Hinc,fvalues(ii),RHBM.r,RHBM.epsilon_r,RHBM.sigma_e,RHBM.rho,tol,RHBM.fN,RHBM.fK);
        % store
        SOL.Jsol(:,:,:,:,1,ii) = Jout;        
        if get(handles.Efield_check, 'Value')
            SOL.Esol(:,:,:,:,1,ii) = Eout;
        end
        if get(handles.Bfield_check, 'Value')
            SOL.Bsol(:,:,:,:,1,ii) = Bout;
        end
        if get(handles.SAR_check, 'Value')
            SOL.Ssol(:,:,:,1,ii) = Sout;
            SOL.Gsar(1,ii) = Gsar;
        end
        if get(handles.Power_check, 'Value')
            SOL.Pabs(1,ii) = Pabs;
        end
    end
    
end


uiwait(msgbox('Scattering Problem Solved','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
uiresume(gcbf);
close(handles.figure1);

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
