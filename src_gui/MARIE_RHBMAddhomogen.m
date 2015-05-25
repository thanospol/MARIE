function varargout = MARIE_RHBMAddhomogen(varargin)
% MARIE_RHBMADDHOMOGEN MATLAB code for MARIE_RHBMAddhomogen.fig
%      MARIE_RHBMADDHOMOGEN, by itself, creates a new MARIE_RHBMADDHOMOGEN or raises the existing
%      singleton*.
%
%      H = MARIE_RHBMADDHOMOGEN returns the handle to a new MARIE_RHBMADDHOMOGEN or the handle to
%      the existing singleton*.
%
%      MARIE_RHBMADDHOMOGEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_RHBMADDHOMOGEN.M with the given input arguments.
%
%      MARIE_RHBMADDHOMOGEN('Property','Value',...) creates a new MARIE_RHBMADDHOMOGEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_RHBMAddhomogen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_RHBMAddhomogen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_RHBMAddhomogen

% Last Modified by GUIDE v2.5 29-Dec-2014 14:45:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_RHBMAddhomogen_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_RHBMAddhomogen_OutputFcn, ...
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


% --- Executes just before MARIE_RHBMAddhomogen is made visible.
function MARIE_RHBMAddhomogen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_RHBMAddhomogen (see VARARGIN)


% Choose default command line output for MARIE_RHBMAddhomogen
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_RHBMAddhomogen wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_RHBMAddhomogen_OutputFcn(hObject, eventdata, handles) 
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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Epsilon_r_Callback(hObject, eventdata, handles)
% hObject    handle to Epsilon_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Epsilon_r as text
%        str2double(get(hObject,'String')) returns contents of Epsilon_r as a double


% --- Executes during object creation, after setting all properties.
function Epsilon_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Epsilon_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Sigma_e_Callback(hObject, eventdata, handles)
% hObject    handle to Sigma_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Sigma_e as text
%        str2double(get(hObject,'String')) returns contents of Sigma_e as a double


% --- Executes during object creation, after setting all properties.
function Sigma_e_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sigma_e (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Density_Callback(hObject, eventdata, handles)
% hObject    handle to Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Density as text
%        str2double(get(hObject,'String')) returns contents of Density as a double


% --- Executes during object creation, after setting all properties.
function Density_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Density (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function Center_X_Cube_Callback(hObject, eventdata, handles)
% hObject    handle to Center_X_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_X_Cube as text
%        str2double(get(hObject,'String')) returns contents of Center_X_Cube as a double


% --- Executes during object creation, after setting all properties.
function Center_X_Cube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_X_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Center_Y_Cube_Callback(hObject, eventdata, handles)
% hObject    handle to Center_Y_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_Y_Cube as text
%        str2double(get(hObject,'String')) returns contents of Center_Y_Cube as a double


% --- Executes during object creation, after setting all properties.
function Center_Y_Cube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_Y_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Center_Z_Cube_Callback(hObject, eventdata, handles)
% hObject    handle to Center_Z_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_Z_Cube as text
%        str2double(get(hObject,'String')) returns contents of Center_Z_Cube as a double


% --- Executes during object creation, after setting all properties.
function Center_Z_Cube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_Z_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Length_X_Cube_Callback(hObject, eventdata, handles)
% hObject    handle to Length_X_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Length_X_Cube as text
%        str2double(get(hObject,'String')) returns contents of Length_X_Cube as a double


% --- Executes during object creation, after setting all properties.
function Length_X_Cube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Length_X_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Length_Y_Cube_Callback(hObject, eventdata, handles)
% hObject    handle to Length_Y_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Length_Y_Cube as text
%        str2double(get(hObject,'String')) returns contents of Length_Y_Cube as a double


% --- Executes during object creation, after setting all properties.
function Length_Y_Cube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Length_Y_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Length_Z_Cube_Callback(hObject, eventdata, handles)
% hObject    handle to Length_Z_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Length_Z_Cube as text
%        str2double(get(hObject,'String')) returns contents of Length_Z_Cube as a double


% --- Executes during object creation, after setting all properties.
function Length_Z_Cube_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Length_Z_Cube (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Add_Cuboid.
function Add_Cuboid_Callback(hObject, eventdata, handles)
% hObject    handle to Add_Cuboid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RHBM;

% get data
Cnt = [ str2num(get(handles.Center_X_Cube, 'String')),...
    str2num(get(handles.Center_Y_Cube, 'String')),...
    str2num(get(handles.Center_Z_Cube, 'String'))];
Len = [ str2num(get(handles.Length_X_Cube, 'String')),...
    str2num(get(handles.Length_Y_Cube, 'String')),...
    str2num(get(handles.Length_Z_Cube, 'String'))];
e_r = str2num(get(handles.Epsilon_r, 'String'));
s_e = str2num(get(handles.Sigma_e, 'String'));
dens = str2num(get(handles.Density, 'String'));


if (sum(find(Len>0)) < 3) || (length(Cnt) < 3)

    % issue warning and stay in function
    f = warndlg('WARNING: wrong cuboid data. Verify all the fields.', 'WRONG DATA');
    waitfor(f);
        
else
    % all ok... create cuboid
    [idx,epsilon_r,sigma_e,rho,~,~] = homogencube(RHBM.r,Cnt,Len,e_r,s_e,dens);
    
    % assign data to RHBM struct... notice it overwrites existing data if
    % overlapping
    RHBM.epsilon_r(idx) = epsilon_r(idx);
    RHBM.sigma_e(idx) = sigma_e(idx);
    RHBM.rho(idx) = rho(idx);
    
    % report and stay
    uiwait(msgbox('Cuboid successfully added to domain','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
    
end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Length_Cyl_Callback(hObject, eventdata, handles)
% hObject    handle to Length_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Length_Cyl as text
%        str2double(get(hObject,'String')) returns contents of Length_Cyl as a double


% --- Executes during object creation, after setting all properties.
function Length_Cyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Length_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ratio_Y_Cyl_Callback(hObject, eventdata, handles)
% hObject    handle to Ratio_Y_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ratio_Y_Cyl as text
%        str2double(get(hObject,'String')) returns contents of Ratio_Y_Cyl as a double


% --- Executes during object creation, after setting all properties.
function Ratio_Y_Cyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ratio_Y_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ratio_X_Cyl_Callback(hObject, eventdata, handles)
% hObject    handle to Ratio_X_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ratio_X_Cyl as text
%        str2double(get(hObject,'String')) returns contents of Ratio_X_Cyl as a double


% --- Executes during object creation, after setting all properties.
function Ratio_X_Cyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ratio_X_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Center_X_Cyl_Callback(hObject, eventdata, handles)
% hObject    handle to Center_X_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_X_Cyl as text
%        str2double(get(hObject,'String')) returns contents of Center_X_Cyl as a double


% --- Executes during object creation, after setting all properties.
function Center_X_Cyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_X_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Center_Y_Cyl_Callback(hObject, eventdata, handles)
% hObject    handle to Center_Y_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_Y_Cyl as text
%        str2double(get(hObject,'String')) returns contents of Center_Y_Cyl as a double


% --- Executes during object creation, after setting all properties.
function Center_Y_Cyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_Y_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Center_Z_Cyl_Callback(hObject, eventdata, handles)
% hObject    handle to Center_Z_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_Z_Cyl as text
%        str2double(get(hObject,'String')) returns contents of Center_Z_Cyl as a double


% --- Executes during object creation, after setting all properties.
function Center_Z_Cyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_Z_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Radius_Cyl_Callback(hObject, eventdata, handles)
% hObject    handle to Radius_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Radius_Cyl as text
%        str2double(get(hObject,'String')) returns contents of Radius_Cyl as a double


% --- Executes during object creation, after setting all properties.
function Radius_Cyl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Radius_Cyl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Add_Cylinder.
function Add_Cylinder_Callback(hObject, eventdata, handles)
% hObject    handle to Add_Cylinder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global RHBM;

% get data
Cnt = [ str2num(get(handles.Center_X_Cyl, 'String')),...
    str2num(get(handles.Center_Y_Cyl, 'String')),...
    str2num(get(handles.Center_Z_Cyl, 'String'))];
Len = str2num(get(handles.Length_Cyl, 'String'));
Rad = str2num(get(handles.Radius_Cyl, 'String'));
Erat = [ str2num(get(handles.Ratio_X_Cyl, 'String')),...
    str2num(get(handles.Ratio_Y_Cyl, 'String'))];
e_r = str2num(get(handles.Epsilon_r, 'String'));
s_e = str2num(get(handles.Sigma_e, 'String'));
dens = str2num(get(handles.Density, 'String'));


if (Len <=0) || (Rad <=0) || (length(Cnt) < 3) || ((sum(find(Erat > 0))) < 2)

    % issue warning and stay in function
    f = warndlg('WARNING: wrong cylinder data. Verify all the fields.', 'WRONG DATA');
    waitfor(f);
        
else
    % all ok... create cylinder
    [idx,epsilon_r,sigma_e,rho,~,~] = homogenellipcylinder(RHBM.r,Cnt,Rad,Len,Erat,e_r,s_e,dens);
    
    % assign data to RHBM struct... notice it overwrites existing data if
    % overlapping
    RHBM.epsilon_r(idx) = epsilon_r(idx);
    RHBM.sigma_e(idx) = sigma_e(idx);
    RHBM.rho(idx) = rho(idx);
    
    % report and stay
    uiwait(msgbox('Cylinder successfully added to domain','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
    
end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function Radius_Ellip_Callback(hObject, eventdata, handles)
% hObject    handle to Radius_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Radius_Ellip as text
%        str2double(get(hObject,'String')) returns contents of Radius_Ellip as a double


% --- Executes during object creation, after setting all properties.
function Radius_Ellip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Radius_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Center_Z_Ellip_Callback(hObject, eventdata, handles)
% hObject    handle to Center_Z_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_Z_Ellip as text
%        str2double(get(hObject,'String')) returns contents of Center_Z_Ellip as a double


% --- Executes during object creation, after setting all properties.
function Center_Z_Ellip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_Z_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Center_Y_Ellip_Callback(hObject, eventdata, handles)
% hObject    handle to Center_Y_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_Y_Ellip as text
%        str2double(get(hObject,'String')) returns contents of Center_Y_Ellip as a double


% --- Executes during object creation, after setting all properties.
function Center_Y_Ellip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_Y_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Center_X_Ellip_Callback(hObject, eventdata, handles)
% hObject    handle to Center_X_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Center_X_Ellip as text
%        str2double(get(hObject,'String')) returns contents of Center_X_Ellip as a double


% --- Executes during object creation, after setting all properties.
function Center_X_Ellip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Center_X_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ratio_X_Ellip_Callback(hObject, eventdata, handles)
% hObject    handle to Ratio_X_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ratio_X_Ellip as text
%        str2double(get(hObject,'String')) returns contents of Ratio_X_Ellip as a double

% --- Executes during object creation, after setting all properties.
function Ratio_X_Ellip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ratio_X_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ratio_Y_Ellip_Callback(hObject, eventdata, handles)
% hObject    handle to Ratio_Y_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ratio_Y_Ellip as text
%        str2double(get(hObject,'String')) returns contents of Ratio_Y_Ellip as a double


% --- Executes during object creation, after setting all properties.
function Ratio_Y_Ellip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ratio_Y_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ratio_Z_Ellip_Callback(hObject, eventdata, handles)
% hObject    handle to Ratio_Z_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ratio_Z_Ellip as text
%        str2double(get(hObject,'String')) returns contents of Ratio_Z_Ellip as a double


% --- Executes during object creation, after setting all properties.
function Ratio_Z_Ellip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ratio_Z_Ellip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Add_Ellipsoid.
function Add_Ellipsoid_Callback(hObject, eventdata, handles)
% hObject    handle to Add_Ellipsoid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global RHBM;

% get data
Cnt = [ str2num(get(handles.Center_X_Ellip, 'String')),...
    str2num(get(handles.Center_Y_Ellip, 'String')),...
    str2num(get(handles.Center_Z_Ellip, 'String'))];
Rad = str2num(get(handles.Radius_Ellip, 'String'));
Erat = [ str2num(get(handles.Ratio_X_Ellip, 'String')),...
    str2num(get(handles.Ratio_Y_Ellip, 'String')),...
    str2num(get(handles.Ratio_Y_Ellip, 'String'))];
e_r = str2num(get(handles.Epsilon_r, 'String'));
s_e = str2num(get(handles.Sigma_e, 'String'));
dens = str2num(get(handles.Density, 'String'));


if (Rad <=0) || (length(Cnt) < 3) || ((sum(find(Erat > 0))) < 3)

    % issue warning and stay in function
    f = warndlg('WARNING: wrong ellipsoid data. Verify all the fields.', 'WRONG DATA');
    waitfor(f);
        
else
    % all ok... create ellipsoid
    [idx,epsilon_r,sigma_e,rho,~,~] = homogenellipse(RHBM.r,Cnt,Rad,Erat,e_r,s_e,dens);
    
    % assign data to RHBM struct... notice it overwrites existing data if
    % overlapping
    RHBM.epsilon_r(idx) = epsilon_r(idx);
    RHBM.sigma_e(idx) = sigma_e(idx);
    RHBM.rho(idx) = rho(idx);   
    RHBM.idxS = find( abs(RHBM.epsilon_r(:) - 1 + 1j*RHBM.sigma_e(:)) > 1e-20);
    
    % report and stay
    uiwait(msgbox('Ellipsoid successfully added to domain','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
    
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcbf);
close(handles.figure1);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Save_Button.
function Save_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RHBM;

% get file name
[status_filename, status_path] = uiputfile('*.vmm', 'Save current status as');

filename = sprintf('%s%s', status_path, status_filename);

% check that return was correct
if (status_filename(1) ~= 0) && (status_path(1) ~= 0)
    fid = fopen(filename, 'w');
    if fid < 1
        f = warndlg('Cannot open selected file', 'Cannot Open File');
        waitfor(f);
    else
        matfilename = status_filename; % change extension to mat
        matfilename(end-3:end)
        matfilename(end-3:end) = '.mat';
        fprintf(fid, '%s\n',RHBM.name);
        fprintf(fid, '%s', matfilename);
        matfilename = sprintf('%s%s', status_path, matfilename);
        r = RHBM.r;
        epsilon_r = RHBM.epsilon_r;
        sigma_e = RHBM.sigma_e;
        rho = RHBM.rho;
        save(matfilename, 'r', 'epsilon_r', 'sigma_e', 'rho', '-v7.3');
        % report and exit
        uiwait(msgbox('Realistic Body Model Saved','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
        uiresume(gcbf);
        close(handles.figure1);
    end
else
    f = warndlg('Incorrect file selected', 'Wrong File');
    waitfor(f);
end
