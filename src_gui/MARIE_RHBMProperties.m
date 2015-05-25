function varargout = MARIE_RHBMProperties(varargin)
% MARIE_RHBMPROPERTIES MATLAB code for MARIE_RHBMProperties.fig
%      MARIE_RHBMPROPERTIES, by itself, creates a new MARIE_RHBMPROPERTIES or raises the existing
%      singleton*.
%
%      H = MARIE_RHBMPROPERTIES returns the handle to a new MARIE_RHBMPROPERTIES or the handle to
%      the existing singleton*.
%
%      MARIE_RHBMPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_RHBMPROPERTIES.M with the given input arguments.
%
%      MARIE_RHBMPROPERTIES('Property','Value',...) creates a new MARIE_RHBMPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_RHBMProperties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_RHBMProperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_RHBMProperties

% Last Modified by GUIDE v2.5 25-Jun-2014 19:54:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_RHBMProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_RHBMProperties_OutputFcn, ...
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


% --- Executes just before MARIE_RHBMProperties is made visible.
function MARIE_RHBMProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_RHBMProperties (see VARARGIN)


% Choose default command line output for MARIE_RHBMProperties
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_RHBMProperties wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_RHBMProperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
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
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


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
if ~isempty(RHBM.name)
    set(hObject,'String',RHBM.name);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Emin_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Emin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Emin as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Emin as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Emin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Emin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.epsilon_r)
    emin = min(unique(RHBM.epsilon_r(:)));
    value = sprintf('%.3g', emin);
    set(hObject,'String',value);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Emax_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Emax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Emax as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Emax as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Emax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Emax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.epsilon_r)
    emax = max(unique(RHBM.epsilon_r(:)));
    value = sprintf('%.3g', emax);
    set(hObject,'String',value);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Smin_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Smin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Smin as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Smin as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Smin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Smin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.sigma_e)
    smin = min(unique(RHBM.sigma_e(:)));
    value = sprintf('%.3g', smin);
    set(hObject,'String',value);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Smax_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Smax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Smax as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Smax as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Smax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Smax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.sigma_e)
    smax = max(unique(RHBM.sigma_e(:)));
    value = sprintf('%.3g', smax);
    set(hObject,'String',value);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Dmin_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Dmin as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Dmin as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Dmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.rho)
    dmin = min(unique(RHBM.rho(:)));
    value = sprintf('%.3g', dmin);
    set(hObject,'String',value);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Dmax_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Dmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Dmax as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Dmax as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Dmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Dmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.rho)
    dmax = max(unique(RHBM.rho(:)));
    value = sprintf('%.3g', dmax);
    set(hObject,'String',value);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Materials_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Materials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Materials as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Materials as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Materials_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Materials (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
mats = [];
if ~isempty(RHBM.epsilon_r)
    mats = RHBM.epsilon_r;
end
if ~isempty(RHBM.sigma_e)
    if ~isempty(mats)
        mats = mats + 1j*RHBM.sigma_e;
    else
        mats = RHBM.sigma_e;
    end
end
if ~isempty(RHBM.rho)
    if ~isempty(mats)
        mats = mats + 1e10*RHBM.rho;
    else
        mats = RHBM.rho;
    end
end
if ~isempty(mats)
    nmats = length(unique(mats(:)));
    value = sprintf('%d', nmats);
    set(hObject,'String',value);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Voxels_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Voxels as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Voxels as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Voxels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
if ~isempty(RHBM.r)
    [L,M,N,~]= size(RHBM.r);
    value = sprintf('%d', 3*L*M*N);
    set(hObject,'String',value);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function RHBM_Unknowns_Callback(hObject, eventdata, handles)
% hObject    handle to RHBM_Unknowns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RHBM_Unknowns as text
%        str2double(get(hObject,'String')) returns contents of RHBM_Unknowns as a double


% --- Executes during object creation, after setting all properties.
function RHBM_Unknowns_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RHBM_Unknowns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global RHBM;
mats = [];
if ~isempty(RHBM.epsilon_r)
    mats = RHBM.epsilon_r - 1;
end
if ~isempty(RHBM.sigma_e)
    if ~isempty(mats)
        mats = mats + 1j*1e6*RHBM.sigma_e;
    else
        mats = 1j*1e6*RHBM.sigma_e;
    end
end
if ~isempty(RHBM.rho)
    if ~isempty(mats)
        mats = mats + 1e6*RHBM.rho;
    else
        mats = 1e6*RHBM.rho;
    end
end
if ~isempty(mats)
    idxS = find(abs(mats(:)) > 1e-12);
    value = sprintf('%d', length(idxS));
    set(hObject,'String',value);
end

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
