function varargout = MARIE_SolveFreq(varargin)
% MARIE_SOLVEFREQ MATLAB code for MARIE_SolveFreq.fig
%      MARIE_SOLVEFREQ, by itself, creates a new MARIE_SOLVEFREQ or raises the existing
%      singleton*.
%
%      H = MARIE_SOLVEFREQ returns the handle to a new MARIE_SOLVEFREQ or the handle to
%      the existing singleton*.
%
%      MARIE_SOLVEFREQ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_SOLVEFREQ.M with the given input arguments.
%
%      MARIE_SOLVEFREQ('Property','Value',...) creates a new MARIE_SOLVEFREQ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_SolveFreq_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_SolveFreq_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_SolveFreq

% Last Modified by GUIDE v2.5 25-Jun-2014 19:36:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_SolveFreq_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_SolveFreq_OutputFcn, ...
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


% --- Executes just before MARIE_SolveFreq is made visible.
function MARIE_SolveFreq_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_SolveFreq (see VARARGIN)

% Choose default command line output for MARIE_SolveFreq
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_SolveFreq wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_SolveFreq_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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
set(hObject,'String',COIL.name);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function FreqMin_Box_Callback(hObject, eventdata, handles)
% hObject    handle to FreqMin_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FreqMin_Box as text
%        str2double(get(hObject,'String')) returns contents of FreqMin_Box as a double


% --- Executes during object creation, after setting all properties.
function FreqMin_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqMin_Box (see GCBO)
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
function FreqMax_Box_Callback(hObject, eventdata, handles)
% hObject    handle to FreqMax_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FreqMax_Box as text
%        str2double(get(hObject,'String')) returns contents of FreqMax_Box as a double


% --- Executes during object creation, after setting all properties.
function FreqMax_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FreqMax_Box (see GCBO)
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
function NumFreqs_Box_Callback(hObject, eventdata, handles)
% hObject    handle to NumFreqs_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumFreqs_Box as text
%        str2double(get(hObject,'String')) returns contents of NumFreqs_Box as a double


% --- Executes during object creation, after setting all properties.
function NumFreqs_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumFreqs_Box (see GCBO)
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
% --- Executes on button press in LinearFreq_Select.
function LinearFreq_Select_Callback(hObject, eventdata, handles)
% hObject    handle to LinearFreq_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LinearFreq_Select

set(handles.LogFreq_Select, 'Value', 0);

% --- Executes on button press in LogFreq_Select.
function LogFreq_Select_Callback(hObject, eventdata, handles)
% hObject    handle to LogFreq_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LogFreq_Select

set(handles.LinearFreq_Select, 'Value', 0);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Solve_Button.
function Solve_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Solve_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get global variables
global COIL;
global SOL;


% get frequency data
numfreqs = round(str2num(get(handles.NumFreqs_Box, 'String')));
freqmin = str2num(get(handles.FreqMin_Box, 'String'));
if (isempty(freqmin))
    errordlg('Wrong Frequency Range -- Please correct', 'Error');
end
if (numfreqs > 1) % frequency sweep
    freqmax = str2num(get(handles.FreqMax_Box, 'String'));
    if (isempty(freqmax)) || (freqmax <= freqmin)
        errordlg('Wrong Frequency Range -- Please correct', 'Error');
    end
    
    if get(handles.LogFreq_Select,'Value') % logscale
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
        
% solve the system
[ZP,Jc] = COIL_Solver(COIL,fvalues);

SOL.Zparam = ZP;
SOL.Jcoil = Jc;
SOL.Jsol = [];
SOL.Esol = [];
SOL.Bsol = [];
SOL.Ssol = [];
SOL.Gsar = [];
SOL.Pabs = [];
SOL.freq = fvalues;

uiwait(msgbox('Frequency Sweep Analysis Done','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
