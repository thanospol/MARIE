function varargout = MARIE_SCOILProperties(varargin)
% MARIE_SCOILPROPERTIES MATLAB code for MARIE_SCOILProperties.fig
%      MARIE_SCOILPROPERTIES, by itself, creates a new MARIE_SCOILPROPERTIES or raises the existing
%      singleton*.
%
%      H = MARIE_SCOILPROPERTIES returns the handle to a new MARIE_SCOILPROPERTIES or the handle to
%      the existing singleton*.
%
%      MARIE_SCOILPROPERTIES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_SCOILPROPERTIES.M with the given input arguments.
%
%      MARIE_SCOILPROPERTIES('Property','Value',...) creates a new MARIE_SCOILPROPERTIES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_SCOILProperties_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_SCOILProperties_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_SCOILProperties

% Last Modified by GUIDE v2.5 23-Jun-2014 14:58:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_SCOILProperties_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_SCOILProperties_OutputFcn, ...
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


% --- Executes just before MARIE_SCOILProperties is made visible.
function MARIE_SCOILProperties_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_SCOILProperties (see VARARGIN)

% Choose default command line output for MARIE_SCOILProperties
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_SCOILProperties wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_SCOILProperties_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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
function Triangles_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Triangles_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Triangles_Box as text
%        str2double(get(hObject,'String')) returns contents of Triangles_Box as a double


% --- Executes during object creation, after setting all properties.
function Triangles_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Triangles_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global COIL;
if ~(strcmp(COIL.name, 'No Selected Model'))
    ntriang = num2str(size(COIL.etod,2));
    set(hObject,'String', ntriang);
end

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
function Nodes_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Nodes_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nodes_Box as text
%        str2double(get(hObject,'String')) returns contents of Nodes_Box as a double


% --- Executes during object creation, after setting all properties.
function Nodes_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nodes_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global COIL;
if ~(strcmp(COIL.name, 'No Selected Model'))
    nnodes = num2str(size(COIL.node,2));
    set(hObject,'String', nnodes);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Unknowns_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Unknowns_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Unknowns_Box as text
%        str2double(get(hObject,'String')) returns contents of Unknowns_Box as a double


% --- Executes during object creation, after setting all properties.
function Unknowns_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Unknowns_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global COIL;
if ~(strcmp(COIL.name, 'No Selected Model'))
    nvars = num2str(max(COIL.index(:)));
    set(hObject,'String', nvars);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Conductivity_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Conductivity_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conductivity_Box as text
%        str2double(get(hObject,'String')) returns contents of Conductivity_Box as a double


% --- Executes during object creation, after setting all properties.
function Conductivity_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conductivity_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global COIL;
if ~(strcmp(COIL.name, 'No Selected Model'))
    if isempty(COIL.Rhocoil)
        rho = 'PEC';
        set(hObject,'String', rho);
    else
        rho = sprintf('%g',1/(COIL.Rhocoil));
        set(hObject,'String', rho);
    end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Edges_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Edges_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edges_Box as text
%        str2double(get(hObject,'String')) returns contents of Edges_Box as a double


% --- Executes during object creation, after setting all properties.
function Edges_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edges_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global COIL;
if ~(strcmp(COIL.name, 'No Selected Model'))
    nedges= num2str(size(COIL.edge,2));
    set(hObject,'String', nedges);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Ports_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Ports_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ports_Box as text
%        str2double(get(hObject,'String')) returns contents of Ports_Box as a double


% --- Executes during object creation, after setting all properties.
function Ports_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ports_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global COIL;
if ~(strcmp(COIL.name, 'No Selected Model'))
    nports = num2str(length(COIL.port));
    set(hObject,'String', nports);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
