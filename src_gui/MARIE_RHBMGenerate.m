function varargout = MARIE_RHBMGenerate(varargin)
% MARIE_RHBMGENERATE MATLAB code for MARIE_RHBMGenerate.fig
%      MARIE_RHBMGENERATE, by itself, creates a new MARIE_RHBMGENERATE or raises the existing
%      singleton*.
%
%      H = MARIE_RHBMGENERATE returns the handle to a new MARIE_RHBMGENERATE or the handle to
%      the existing singleton*.
%
%      MARIE_RHBMGENERATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_RHBMGENERATE.M with the given input arguments.
%
%      MARIE_RHBMGENERATE('Property','Value',...) creates a new MARIE_RHBMGENERATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_RHBMGenerate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_RHBMGenerate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_RHBMGenerate

% Last Modified by GUIDE v2.5 29-Dec-2014 11:18:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_RHBMGenerate_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_RHBMGenerate_OutputFcn, ...
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


% --- Executes just before MARIE_RHBMGenerate is made visible.
function MARIE_RHBMGenerate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_RHBMGenerate (see VARARGIN)


% Choose default command line output for MARIE_RHBMGenerate
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_RHBMGenerate wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_RHBMGenerate_OutputFcn(hObject, eventdata, handles) 
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


% create the r and main data of the new domain
% limits of domain
applydomain = 0;
if ~isempty(xlimits)
    if (nnz(xlimits))
        xmin = min(xlimits);
        xmax = max(xlimits);
        if (xmin ~= xmax)
            applydomain = applydomain + 1;
        end
    end
end

if ~isempty(ylimits)
    if (nnz(ylimits))
        ymin = min(ylimits);
        ymax = max(ylimits);
        if (ymin ~= ymax)
            applydomain = applydomain + 1;
        end
    end
end

if ~isempty(zlimits)
    if (nnz(zlimits))
        zmin = min(zlimits);
        zmax = max(zlimits);
        if (zmin ~= zmax)
            applydomain = applydomain + 1;
        end
    end
end

% define discretization
if ~isempty(resolution)
    % force resolution
    dx = resolution;
    applydomain = applydomain + 1;
end

if (applydomain < 4) % error: not enough data
    
    % issue warning and stay in function
    f = warndlg('WARNING: wrong data. Verify all the fields.', 'WRONG DATA');
    waitfor(f);
    
else
    % number of voxels in each direction
    nX = ceil((xmax - xmin)/dx);
    nY = ceil((ymax - ymin)/dx);
    nZ = ceil((zmax - zmin)/dx);
    
    % make sure we cover the whole domain (extra voxel in case non exact)
    x = xmin:dx:xmin+nX*dx;
    y = ymin:dx:ymin+nY*dx;
    z = zmin:dx:zmin+nZ*dx;
    
    % generate the 3D grid
    r = grid3d(x,y,z);
    [L,M,N,~] = size(r);
    
    % generate initial (empty) data
    epsilon_r = ones(L,M,N);
    sigma_e = zeros(L,M,N);
    rho = zeros(L,M,N);
    idxS = [];
    
    % name
    name = get(handles.RHBM_ModelName, 'String');
    
    % Assign data to structure
    RHBM = struct('name', name,...
        'r',  r,...
        'epsilon_r', epsilon_r, ...
        'sigma_e', sigma_e, ...
        'rho', rho, ...
        'idxS', idxS, ...
        'freqfN', [], ...
        'freqfK', [], ...
        'fN', [], ...
        'fK', [], ...
        'freqM', [], ...
        'Dcoord', [], ...
        'P', [], ...
        'Um', [], ...
        'Sm', [], ...
        'Vm', [], ...
        'X', [], ...
        'M', []);
    
    % call the function to add homogeneous data
    MARIE_RHBMAddhomogen;
    uiwait(gcf);
    
    % message and close
    % uiwait(msgbox('RHBM generated','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
    uiresume(gcbf);
    close(handles.figure1);
    
end


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




% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcbf);
close(handles.figure1);
