function varargout = MARIE_ExportSnP(varargin)
% MARIE_EXPORTSNP MATLAB code for MARIE_ExportSnP.fig
%      MARIE_EXPORTSNP, by itself, creates a new MARIE_EXPORTSNP or raises the existing
%      singleton*.
%
%      H = MARIE_EXPORTSNP returns the handle to a new MARIE_EXPORTSNP or the handle to
%      the existing singleton*.
%
%      MARIE_EXPORTSNP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_EXPORTSNP.M with the given input arguments.
%
%      MARIE_EXPORTSNP('Property','Value',...) creates a new MARIE_EXPORTSNP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_ExportSnP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_ExportSnP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_ExportSnP

% Last Modified by GUIDE v2.5 10-Nov-2014 17:12:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_ExportSnP_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_ExportSnP_OutputFcn, ...
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


% --- Executes just before MARIE_ExportSnP is made visible.
function MARIE_ExportSnP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_ExportSnP (see VARARGIN)

% Choose default command line output for MARIE_ExportSnP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_ExportSnP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_ExportSnP_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Zparam_button.
function Zparam_button_Callback(hObject, eventdata, handles)
% hObject    handle to Zparam_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Zparam_button

set(handles.Sparam_button, 'Value', 0);
set(handles.Yparam_button, 'Value', 0);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Yparam_button.
function Yparam_button_Callback(hObject, eventdata, handles)
% hObject    handle to Yparam_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Yparam_button

set(handles.Zparam_button, 'Value', 0);
set(handles.Sparam_button, 'Value', 0);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Sparam_button.
function Sparam_button_Callback(hObject, eventdata, handles)
% hObject    handle to Sparam_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sparam_button

set(handles.Zparam_button, 'Value', 0);
set(handles.Yparam_button, 'Value', 0);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Refload_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Refload_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Refload_Edit as text
%        str2double(get(hObject,'String')) returns contents of Refload_Edit as a double


% --- Executes during object creation, after setting all properties.
function Refload_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Refload_Edit (see GCBO)
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
% --- Executes on button press in Export_button.
function Export_button_Callback(hObject, eventdata, handles)
% hObject    handle to Export_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SOL;

Nports = size(SOL.Zparam,1);

% type of parameters
Type = [];
if get(handles.Yparam_button,'Value') % Y parameter
    Type = 'Y';
end
if get(handles.Zparam_button,'Value') % Z parameter
    Type = 'Z';
end
if get(handles.Sparam_button,'Value') % S parameter
    Type = 'S';
end

% reference load
Zref = [];
Zref = str2num(get(handles.Refload_Edit, 'String'));
if (Zref <= 0)
    Zref = [];
end


% get file name
[snp_filename, snp_path] = uiputfile(sprintf('*.s%dp',Nports), 'Export Port Parameters to');

filename = sprintf('%s%s', snp_path, snp_filename);

% check that return was correct
if (snp_filename(1) ~= 0) && (snp_path(1) ~= 0)
    Export_SnP(filename,SOL.Zparam,SOL.freq,Type,Zref);
else
    f = warndlg('Incorrect file selected', 'Wrong File');
    waitfor(f);
end
    
uiwait(msgbox('Port Parameters saved to Touchstone format file','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
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
