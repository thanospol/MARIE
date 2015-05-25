function varargout = MARIE_SolveMR(varargin)
% MARIE_SOLVEMR MATLAB code for MARIE_SolveMR.fig
%      MARIE_SOLVEMR, by itself, creates a new MARIE_SOLVEMR or raises the existing
%      singleton*.
%
%      H = MARIE_SOLVEMR returns the handle to a new MARIE_SOLVEMR or the handle to
%      the existing singleton*.
%
%      MARIE_SOLVEMR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_SOLVEMR.M with the given input arguments.
%
%      MARIE_SOLVEMR('Property','Value',...) creates a new MARIE_SOLVEMR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_SolveMR_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_SolveMR_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_SolveMR

% Last Modified by GUIDE v2.5 14-May-2015 15:35:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_SolveMR_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_SolveMR_OutputFcn, ...
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


% --- Executes just before MARIE_SolveMR is made visible.
function MARIE_SolveMR_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_SolveMR (see VARARGIN)

% Choose default command line output for MARIE_SolveMR
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_SolveMR wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_SolveMR_OutputFcn(hObject, eventdata, handles) 
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
% --- Executes on button press in Solve_Button.
function Solve_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Solve_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get global variables
global RHBM;
global COIL;
global SOL;


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
        
% get tolerance
tol = str2num(get(handles.Tolerance_Edit, 'String'));

% set the values to compute
flags = zeros(9,1);
flags(2) = 1; % compute the coil currents by default

if strcmp(RHBM.name,'No Model Selected');

    validcoil = 1;

    switch COIL.type
        case 'S'
            
            % get positions of the coil elements
            xmin = min(COIL.node(:,1));
            xmax = max(COIL.node(:,1));
            ymin = min(COIL.node(:,1));
            ymax = max(COIL.node(:,1));
            zmin = min(COIL.node(:,1));
            zmax = max(COIL.node(:,1));
            
        case 'W'
            
            % get positions of the coil elements
            xmin = min(COIL.Pcoil(:,1));
            xmax = max(COIL.Pcoil(:,1));
            ymin = min(COIL.Pcoil(:,2));
            ymax = max(COIL.Pcoil(:,2));
            zmin = min(COIL.Pcoil(:,3));
            zmax = max(COIL.Pcoil(:,3));

        otherwise
            
            validcoil = 0;
            
    end

    if validcoil
        % create a new grid
        dx = 0.005; % 5mm grid
        % number of voxels in each direction
        nX = ceil((xmax - xmin)/dx)+1;
        nY = ceil((ymax - ymin)/dx)+1;
        nZ = ceil((zmax - zmin)/dx)+1;
        
        % make sure we cover the whole domain (extra voxel in case non exact)
        xnew = xmin-dx:dx:xmin+nX*dx;
        ynew = ymin-dx:dx:ymin+nY*dx;
        znew = zmin-dx:dx:zmin+nZ*dx;
        
        % generate the 3D grid and assign free-space to the grid
        RHBM.r = grid3d(xnew,ynew,znew);
        [L,M,N,~] = size(RHBM.r);
        RHBM.epsilon_r = ones(L,M,N);
        RHBM.sigma_e = zeros(L,M,N);
        RHBM.rho = zeros(L,M,N);
        RHBM.name = 'Free Space';
        
    else
        flags(2) = 0;        
    end

else
    % we have RHBM
    flags(3) = 1;

end

if get(handles.Sparam_Select, 'Value')
    flags(1) = 1;
end
if get(handles.EField_Select, 'Value')
    flags(5) = 1;
    flags(8) = 1;
end
if get(handles.BField_Select, 'Value')
    flags(6) = 1;
end
if get(handles.SAR_Select, 'Value')
    flags(4) = 1;
    flags(7) = 1;
end

if get(handles.InAir_Select, 'Value')
    flags(9) = 1;
end

% solve the system
[ZP,Jc,Jb,Sb,Eb,Bb,Gsar,Pabs] = MR_Solver(RHBM,COIL,fvalues,tol,flags);

SOL.Zparam = ZP;
SOL.Jcoil = Jc;
SOL.Jsol = Jb;
SOL.Esol = Eb;
SOL.Bsol = Bb;
SOL.Ssol = Sb;
SOL.Gsar = Gsar;
SOL.Pabs = Pabs;
SOL.freq = fvalues;


uiwait(msgbox('Full-Wave Analysis Done','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in EField_Select.
function EField_Select_Callback(hObject, eventdata, handles)
% hObject    handle to EField_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of EField_Select


% --- Executes on button press in BField_Select.
function BField_Select_Callback(hObject, eventdata, handles)
% hObject    handle to BField_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BField_Select


% --- Executes on button press in InAir_Select.
function InAir_Select_Callback(hObject, eventdata, handles)
% hObject    handle to InAir_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of InAir_Select


% --- Executes on button press in SAR_Select.
function SAR_Select_Callback(hObject, eventdata, handles)
% hObject    handle to SAR_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SAR_Select


% --- Executes on button press in Sparam_Select.
function Sparam_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Sparam_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sparam_Select

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

function Tolerance_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Tolerance_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tolerance_Edit as text
%        str2double(get(hObject,'String')) returns contents of Tolerance_Edit as a double


% --- Executes during object creation, after setting all properties.
function Tolerance_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tolerance_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
