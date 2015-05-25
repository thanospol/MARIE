function varargout = MARIE_SCOILGenerate(varargin)
% MARIE_SCOILGENERATE MATLAB code for MARIE_SCOILGenerate.fig
%      MARIE_SCOILGENERATE, by itself, creates a new MARIE_SCOILGENERATE or raises the existing
%      singleton*.
%
%      H = MARIE_SCOILGENERATE returns the handle to a new MARIE_SCOILGENERATE or the handle to
%      the existing singleton*.
%
%      MARIE_SCOILGENERATE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_SCOILGENERATE.M with the given input arguments.
%
%      MARIE_SCOILGENERATE('Property','Value',...) creates a new MARIE_SCOILGENERATE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_SCOILGenerate_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_SCOILGenerate_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_SCOILGenerate

% Last Modified by GUIDE v2.5 01-Oct-2014 16:07:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_SCOILGenerate_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_SCOILGenerate_OutputFcn, ...
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


% --- Executes just before MARIE_SCOILGenerate is made visible.
function MARIE_SCOILGenerate_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_SCOILGenerate (see VARARGIN)


% get initial variables
handles.coilnumber = varargin{1};
handles.folder = varargin{2}; % current DATA folder

% Choose default command line output for MARIE_SCOILGenerate
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_SCOILGenerate wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_SCOILGenerate_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Modelname_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Modelname_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Modelname_Edit as text
%        str2double(get(hObject,'String')) returns contents of Modelname_Edit as a double


% --- Executes during object creation, after setting all properties.
function Modelname_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Modelname_Edit (see GCBO)
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
function Xposition_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Xposition_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xposition_Edit as text
%        str2double(get(hObject,'String')) returns contents of Xposition_Edit as a double


% --- Executes during object creation, after setting all properties.
function Xposition_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xposition_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Yposition_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Yposition_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Yposition_Edit as text
%        str2double(get(hObject,'String')) returns contents of Yposition_Edit as a double


% --- Executes during object creation, after setting all properties.
function Yposition_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Yposition_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Zposition_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Zposition_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Zposition_Edit as text
%        str2double(get(hObject,'String')) returns contents of Zposition_Edit as a double


% --- Executes during object creation, after setting all properties.
function Zposition_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zposition_Edit (see GCBO)
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
function Coiltrace_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Coiltrace_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Coiltrace_Edit as text
%        str2double(get(hObject,'String')) returns contents of Coiltrace_Edit as a double


% --- Executes during object creation, after setting all properties.
function Coiltrace_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Coiltrace_Edit (see GCBO)
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
function Portgap_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Portgap_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Portgap_Edit as text
%        str2double(get(hObject,'String')) returns contents of Portgap_Edit as a double


% --- Executes during object creation, after setting all properties.
function Portgap_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Portgap_Edit (see GCBO)
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
function Coilradius_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Coilradius_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Coilradius_Edit as text
%        str2double(get(hObject,'String')) returns contents of Coilradius_Edit as a double


% --- Executes during object creation, after setting all properties.
function Coilradius_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Coilradius_Edit (see GCBO)
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
function Shield_Length_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Shield_Length_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Shield_Length_Edit as text
%        str2double(get(hObject,'String')) returns contents of Shield_Length_Edit as a double


% --- Executes during object creation, after setting all properties.
function Shield_Length_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Shield_Length_Edit (see GCBO)
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
function Shield_Radius_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Shield_Radius_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Shield_Radius_Edit as text
%        str2double(get(hObject,'String')) returns contents of Shield_Radius_Edit as a double


% --- Executes during object creation, after setting all properties.
function Shield_Radius_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Shield_Radius_Edit (see GCBO)
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
function Conformal_Numbercoils_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_Numbercoils_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conformal_Numbercoils_Edit as text
%        str2double(get(hObject,'String')) returns contents of Conformal_Numbercoils_Edit as a double


% --- Executes during object creation, after setting all properties.
function Conformal_Numbercoils_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conformal_Numbercoils_Edit (see GCBO)
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
function Conformal_Coillength_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_Coillength_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conformal_Coillength_Edit as text
%        str2double(get(hObject,'String')) returns contents of Conformal_Coillength_Edit as a double


% --- Executes during object creation, after setting all properties.
function Conformal_Coillength_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conformal_Coillength_Edit (see GCBO)
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
function Conformal_Coilspan_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_Coilspan_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conformal_Coilspan_Edit as text
%        str2double(get(hObject,'String')) returns contents of Conformal_Coilspan_Edit as a double


% --- Executes during object creation, after setting all properties.
function Conformal_Coilspan_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conformal_Coilspan_Edit (see GCBO)
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
function Conformal_Iniangle_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_Iniangle_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conformal_Iniangle_Edit as text
%        str2double(get(hObject,'String')) returns contents of Conformal_Iniangle_Edit as a double


% --- Executes during object creation, after setting all properties.
function Conformal_Iniangle_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conformal_Iniangle_Edit (see GCBO)
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
function Planar_Numbercoils_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_Numbercoils_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Planar_Numbercoils_Edit as text
%        str2double(get(hObject,'String')) returns contents of Planar_Numbercoils_Edit as a double


% --- Executes during object creation, after setting all properties.
function Planar_Numbercoils_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Planar_Numbercoils_Edit (see GCBO)
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
function Planar_Coillength_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_Coillength_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Planar_Coillength_Edit as text
%        str2double(get(hObject,'String')) returns contents of Planar_Coillength_Edit as a double


% --- Executes during object creation, after setting all properties.
function Planar_Coillength_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Planar_Coillength_Edit (see GCBO)
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
function Planar_Coilspan_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_Coilspan_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Planar_Coilspan_Edit as text
%        str2double(get(hObject,'String')) returns contents of Planar_Coilspan_Edit as a double


% --- Executes during object creation, after setting all properties.
function Planar_Coilspan_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Planar_Coilspan_Edit (see GCBO)
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
function Planar_Iniangle_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_Iniangle_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Planar_Iniangle_Edit as text
%        str2double(get(hObject,'String')) returns contents of Planar_Iniangle_Edit as a double


% --- Executes during object creation, after setting all properties.
function Planar_Iniangle_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Planar_Iniangle_Edit (see GCBO)
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
function Bird_Numberlegs_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Bird_Numberlegs_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bird_Numberlegs_Edit as text
%        str2double(get(hObject,'String')) returns contents of Bird_Numberlegs_Edit as a double


% --- Executes during object creation, after setting all properties.
function Bird_Numberlegs_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bird_Numberlegs_Edit (see GCBO)
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
function Bird_Coillength_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Bird_Coillength_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bird_Coillength_Edit as text
%        str2double(get(hObject,'String')) returns contents of Bird_Coillength_Edit as a double


% --- Executes during object creation, after setting all properties.
function Bird_Coillength_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bird_Coillength_Edit (see GCBO)
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
function Bird_Iniangle_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Bird_Iniangle_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bird_Iniangle_Edit as text
%        str2double(get(hObject,'String')) returns contents of Bird_Iniangle_Edit as a double


% --- Executes during object creation, after setting all properties.
function Bird_Iniangle_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bird_Iniangle_Edit (see GCBO)
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
function Planar_Topport_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_Topport_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Planar_Topport_Edit as text
%        str2double(get(hObject,'String')) returns contents of Planar_Topport_Edit as a double


% --- Executes during object creation, after setting all properties.
function Planar_Topport_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Planar_Topport_Edit (see GCBO)
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
function Planar_Rightport_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_Rightport_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Planar_Rightport_Edit as text
%        str2double(get(hObject,'String')) returns contents of Planar_Rightport_Edit as a double


% --- Executes during object creation, after setting all properties.
function Planar_Rightport_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Planar_Rightport_Edit (see GCBO)
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
function Planar_Bottomport_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_Bottomport_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Planar_Bottomport_Edit as text
%        str2double(get(hObject,'String')) returns contents of Planar_Bottomport_Edit as a double


% --- Executes during object creation, after setting all properties.
function Planar_Bottomport_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Planar_Bottomport_Edit (see GCBO)
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
function Planar_Leftport_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_Leftport_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Planar_Leftport_Edit as text
%        str2double(get(hObject,'String')) returns contents of Planar_Leftport_Edit as a double


% --- Executes during object creation, after setting all properties.
function Planar_Leftport_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Planar_Leftport_Edit (see GCBO)
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
function Conformal_Topport_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_Topport_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conformal_Topport_Edit as text
%        str2double(get(hObject,'String')) returns contents of Conformal_Topport_Edit as a double


% --- Executes during object creation, after setting all properties.
function Conformal_Topport_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conformal_Topport_Edit (see GCBO)
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
function Conformal_Rightport_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_Rightport_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conformal_Rightport_Edit as text
%        str2double(get(hObject,'String')) returns contents of Conformal_Rightport_Edit as a double


% --- Executes during object creation, after setting all properties.
function Conformal_Rightport_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conformal_Rightport_Edit (see GCBO)
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
function Conformal_Bottomport_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_Bottomport_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conformal_Bottomport_Edit as text
%        str2double(get(hObject,'String')) returns contents of Conformal_Bottomport_Edit as a double


% --- Executes during object creation, after setting all properties.
function Conformal_Bottomport_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conformal_Bottomport_Edit (see GCBO)
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
function Conformal_Leftport_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_Leftport_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Conformal_Leftport_Edit as text
%        str2double(get(hObject,'String')) returns contents of Conformal_Leftport_Edit as a double


% --- Executes during object creation, after setting all properties.
function Conformal_Leftport_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Conformal_Leftport_Edit (see GCBO)
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
% --- Executes on button press in Conformal_button.
function Conformal_button_Callback(hObject, eventdata, handles)
% hObject    handle to Conformal_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Conformal_button

set(handles.Planar_button, 'Value', 0);
set(handles.Bird_button, 'Value', 0);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Planar_button.
function Planar_button_Callback(hObject, eventdata, handles)
% hObject    handle to Planar_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Planar_button

set(handles.Conformal_button, 'Value', 0);
set(handles.Bird_button, 'Value', 0);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Bird_button.
function Bird_button_Callback(hObject, eventdata, handles)
% hObject    handle to Bird_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Bird_button

set(handles.Planar_button, 'Value', 0);
set(handles.Conformal_button, 'Value', 0);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Shield_button.
function Shield_button_Callback(hObject, eventdata, handles)
% hObject    handle to Shield_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Shield_button


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Generate_Button.
function Generate_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Generate_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get name
coilname = get(handles.Modelname_Edit, 'String');
coilpath = handles.folder;
filename = strcat(coilpath, coilname);

% get position data
Xpos = str2num(get(handles.Xposition_Edit, 'String'));
Ypos = str2num(get(handles.Yposition_Edit, 'String'));
Zpos = str2num(get(handles.Zposition_Edit, 'String'));


% get coil info 
Rad = str2num(get(handles.Coilradius_Edit, 'String'));
Wid = str2num(get(handles.Coiltrace_Edit, 'String'));
Gap = str2num(get(handles.Portgap_Edit, 'String'));


% get shield info
shieldflag = get(handles.Shield_button, 'Value');
if shieldflag
    ShieldRad = str2num(get(handles.Shield_Radius_Edit, 'String'));
    ShieldLen = str2num(get(handles.Shield_Length_Edit, 'String'));
else
    ShieldRad = [];
    ShieldLen = [];
end
    

% get type of coil
type = find([get(handles.Conformal_button, 'Value') get(handles.Planar_button, 'Value') get(handles.Bird_button, 'Value')]);

switch type

case 1
    % conformal coil case

    % get number of coils
    Ncoils  = str2num(get(handles.Conformal_Numbercoils_Edit, 'String'));
    % get length of coil
    Len  = str2num(get(handles.Conformal_Coillength_Edit, 'String'));
    % get span of each coil
    Asp  = str2num(get(handles.Conformal_Coilspan_Edit, 'String'));
    % get initial rotation angle
    Iniangle = str2num(get(handles.Conformal_Iniangle_Edit, 'String'));
    % get top ports
    Tports = str2num(get(handles.Conformal_Topport_Edit, 'String'));
    % get bottom ports
    Bports = str2num(get(handles.Conformal_Bottomport_Edit, 'String'));
    % get left ports
    Lports = str2num(get(handles.Conformal_Leftport_Edit, 'String'));
    % get right ports
    Rports = str2num(get(handles.Conformal_Rightport_Edit, 'String'));
    
    % call the coil generator
    Conformal_SCOIL_Gen(filename,Ncoils,Tports,Rports,Bports,Lports,Rad,Xpos,Ypos,Zpos,Len,Asp,Wid,Gap,Iniangle,ShieldRad,ShieldLen);
    
    % generate the file information in the scm file
    scoilfile = sprintf('%s.smm',filename);
    fid = fopen(scoilfile,'w');
    fprintf(fid,'Conformal %d coil array, R = %f m, Length = %f m \n', Ncoils, Rad, Len);
    Rho = 1.2579e-08; thickness = 0.0005;
    fprintf(fid,'%f %f\n', Rho, thickness);
    fprintf(fid,'%s.msh', coilname);
    fclose(fid);

case 2
    % planar coil case

    % get number of coils
    Ncoils  = str2num(get(handles.Planar_Numbercoils_Edit, 'String'));
    % get length of coil
    Len  = str2num(get(handles.Planar_Coillength_Edit, 'String'));
    % get span of each coil
    Ssp  = str2num(get(handles.Planar_Coilspan_Edit, 'String'));
    % get initial rotation angle
    Iniangle = str2num(get(handles.Planar_Iniangle_Edit, 'String'));
    % get top ports
    Tports = str2num(get(handles.Planar_Topport_Edit, 'String'));
    % get bottom ports
    Bports = str2num(get(handles.Planar_Bottomport_Edit, 'String'));
    % get left ports
    Lports = str2num(get(handles.Planar_Leftport_Edit, 'String'));
    % get right ports
    Rports = str2num(get(handles.Planar_Rightport_Edit, 'String'));
    
    % call the coil generator
    Planar_SCOIL_Gen(filename,Ncoils,Tports,Rports,Bports,Lports,Rad,Xpos,Ypos,Zpos,Len,Ssp,Wid,Gap,Iniangle,ShieldRad,ShieldLen);

    % generate the file information in the scm file
    scoilfile = sprintf('%s.smm',filename);
    fid = fopen(scoilfile,'w');
    fprintf(fid,'Planar %d coil array, R = %f m, Length = %f m \n', Ncoils, Rad, Len);
    Rho = 1.2579e-08; thickness = 0.0005;
    fprintf(fid,'%f %f\n', Rho, thickness);
    fprintf(fid,'%s.msh', coilname);
    fclose(fid);   
    

case 3
    % birdcage

    % get number of coils
    Nleg  = str2num(get(handles.Bird_Numberlegs_Edit, 'String'));
    % get length of coil
    Len  = str2num(get(handles.Bird_Coillength_Edit, 'String'));
    % get initial rotation angle
    Iniangle = str2num(get(handles.Bird_Iniangle_Edit, 'String'));

    % call the generator
    Birdcage_SCOIL_Gen(filename,Nleg,Rad,Xpos,Ypos,Zpos,Len,Wid,Gap,Iniangle,ShieldRad,ShieldLen);
    
    % generate the file information in the scm file
    scoilfile = sprintf('%s.smm',filename);
    fid = fopen(scoilfile,'w');
    fprintf(fid,'Birdcage coil, R = %f m, Length = %f m \n', Rad, Len);
    Rho = 1.2579e-08; thickness = 0.0005;
    fprintf(fid,'%f %f\n', Rho, thickness);
    fprintf(fid,'%s.msh', coilname);
    fclose(fid);    
    
    
otherwise
    % No coil

end

meshfile = sprintf('%s.msh', coilname);
if exist(meshfile, 'file')
    % the mesh file was generated correctly, all OK
    uiwait(msgbox('Coil Generated Correctly.','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
    uiresume(gcbf);
else
    % the mesh file was not generated correctly
    f = warndlg('WARNING: .msh file not generated, coil invalid. Make sure to generate the  mesh in gmsh (check documentation for details).', 'MESH INCOMPLETE');
    waitfor(f);
end
close(handles.figure1);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Closebutton.
function Closebutton_Callback(hObject, eventdata, handles)
% hObject    handle to Closebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcbf);
close(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
