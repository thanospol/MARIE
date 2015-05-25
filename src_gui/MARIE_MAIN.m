function varargout = MARIE_MAIN(varargin)
% MARIE_MAIN MATLAB code for MARIE_MAIN.fig
%      MARIE_MAIN, by itself, creates a new MARIE_MAIN or raises the existing
%      singleton*.
%
%      H = MARIE_MAIN returns the handle to a new MARIE_MAIN or the handle to
%      the existing singleton*.
%
%      MARIE_MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_MAIN.M with the given input arguments.
%
%      MARIE_MAIN('Property','Value',...) creates a new MARIE_MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_MAIN_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_MAIN_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_MAIN

% Last Modified by GUIDE v2.5 10-Nov-2014 16:51:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_MAIN_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_MAIN_OutputFcn, ...
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

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes just before MARIE_MAIN is made visible.
function MARIE_MAIN_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_MAIN (see VARARGIN)

% Choose default command line output for MARIE_MAIN
handles.output = hObject;


handles.Ncoils = 0;
handles.Nbodies = 0;
if ispc
   % define the path to the RHBM
    handles.path = strcat(pwd, '\data\');
else
    % define the path to the RHBM
    handles.path = strcat(pwd, '/data/');
end

% initialize RHBM global variable
global RHBM;
RHBM.name = 'No Model Selected';
RHBM.r = [];
RHBM.epsilon_r = [];
RHBM.sigma_e = [];
RHBM.rho = [];

% initialize COIL global variable
global COIL;
COIL.name = 'No Model Selected';
COIL.type = [];

% initialize solution
global SOL;
SOL.Jsol = [];



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_MAIN wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = MARIE_MAIN_OutputFcn(hObject, eventdata, handles) 
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
% --- Executes during object creation, after setting all properties.
function Logo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Logo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Logo
FIG = gcbf;
figure(FIG);

if ispc
    % create the LOGO
    logoim = imread('MarieFulllogo_blueorange.jpg');
    image(logoim);
    axis image;
    set(gca,'XTick', []);
    set(gca,'YTick', []);
    
else
    % create the LOGO
    logoim = imread('MarieFulllogo_blueorange.jpg');
    image(logoim);
    axis image;
    set(gca,'XTick', []);
    set(gca,'YTick', []);

end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function Geometry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Geometry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate Geometry
view(3);
xlabel(gca,'x (m)','Fontsize', 14);
ylabel(gca,'y (m)','Fontsize', 14);
zlabel(gca,'z (m)','Fontsize', 14);
set(gca, 'DataAspectRatio', [1,1,1]);
axis equal;
grid on;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Exit.
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

q = questdlg('Exit MARIE?','Exit and Close','Yes','No','No');
if strcmp(q,'Yes' )
    close all hidden;
    clear all;
    clc;
    fprintf(1, '\n\n');
    fprintf(1, '-----------------------------------------------------------\n');
    fprintf(1, '  MARIE - Magnetic Resonance Integral Equation suite\n');
    fprintf(1, '          Jorge Fernandez Villena   -- jvillena@mit.edu\n');
    fprintf(1, '          Athanasios G. Polimeridis -- thanos_p@mit.edu\n');
    fprintf(1, '          Copyright © 2014\n');
    fprintf(1, '          RLE Computational Prototyping Group, MIT\n\n');
    fprintf(1, '          This software is free and open source\n');
    fprintf(1, '          Distributed under the GNU-GPLv3 terms\n');
    fprintf(1, '          For details see MARIE_license.txt\n\n');
    fprintf(1, '          Thanks for using MARIE!\n');
    fprintf(1, '-----------------------------------------------------------\n\n\n');
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in License.
function License_Callback(hObject, eventdata, handles)
% hObject    handle to License (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_License; % call the license display
uiwait(gcf);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Refresh.
function Refresh_Callback(hObject, eventdata, handles)
% hObject    handle to Refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RHBM;
global COIL;

% call the function to update the Geometry Plot
str = get(handles.Plot_ONOFF_Button,'String');    % get the current pushbutton string 
Plot_on = strcmp(str,'Plot ON'); 
Azimut = get(handles.Azimut_Slider,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
update_Geometry(handles.Geometry,Azimut,Elevation,Plot_on);

% update the entries of the models
if ~isempty(RHBM.name)
    set(handles.ModelName_RHBM,'String',RHBM.name);
end
if ~isempty(COIL.name)
    set(handles.ModelName_COIL,'String',COIL.name);
end
% Commit the new struct element to appdata
guidata(hObject, handles); 

clc;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Clear_Button.
function Clear_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Clear_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear RHBM;
clear COIL;
clear SOL;
clear PULSE;

% initialize RHBM global variable
global RHBM;
RHBM.name = 'No Model Selected';
RHBM.r = [];
RHBM.epsilon_r = [];
RHBM.sigma_e = [];
RHBM.rho = [];
RHBM.idxS = [];

% initialize COIL global variable
global COIL;
COIL.name = 'No Model Selected';
COIL.type = 'N';

% initialize solution
global SOL;
SOL.Jsol = [];
SOL.Esol = [];
SOL.Bsol = [];
SOL.Ssol = [];
SOL.Gsar = [];
SOL.Pabs = [];
SOL.Sparam = [];
SOL.freq = [];

global PULSE;
PULSE.Jcoil = [];
PULSE.Jsol = [];
PULSE.Esol = [];
PULSE.Bsol = [];
PULSE.Ssol = [];
PULSE.Gsar = [];
PULSE.Pabs = [];
PULSE.Ipulse = [];
PULSE.Vpulse = [];
PULSE.freq = [];


% call the function to update the Geometry Plot
Azimut = get(handles.Azimut_Slider,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
update_Geometry(handles.Geometry,Azimut,Elevation,0);

% update the entries of the models
if ~isempty(RHBM.name)
    set(handles.ModelName_RHBM,'String',RHBM.name);
end
if ~isempty(COIL.name)
    set(handles.ModelName_COIL,'String',COIL.name);
end
% Commit the new struct element to appdata
guidata(hObject, handles); 

clc;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Simulate.
function Simulate_Callback(hObject, eventdata, handles)
% hObject    handle to Simulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_Analyze;
uiwait(gcf);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Visualize.
function Visualize_Callback(hObject, eventdata, handles)
% hObject    handle to Visualize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_Visualize;
uiwait(gcf);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Properties_COIL.
function Properties_COIL_Callback(hObject, eventdata, handles)
% hObject    handle to Properties_COIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global COIL;
if ~isempty(COIL.type)
    switch COIL.type        
        case 'W'
            MARIE_WCOILProperties;
        case 'S'
            MARIE_SCOILProperties;
    end
    uiwait(gcf);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Modify_COIL.
function Modify_COIL_Callback(hObject, eventdata, handles)
% hObject    handle to Modify_COIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_COILModify(handles.Geometry); % call the modify GUI
uiwait(gcf);

% call the function to update the Geometry Plot
str = get(handles.Plot_ONOFF_Button,'String');    % get the current pushbutton string 
Plot_on = strcmp(str,'Plot ON'); 
Azimut = get(handles.Azimut_Slider,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
update_Geometry(handles.Geometry,Azimut,Elevation,Plot_on);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Select_COIL.
function Select_COIL_Callback(hObject, eventdata, handles)
% hObject    handle to Select_COIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



if ispc
   % define the path to the coil
    coilpath = strcat(handles.path, 'coils\');
else
    % define the path to the coil
    coilpath = strcat(handles.path, 'coils/');
end

MARIE_COILLoad(handles.Ncoils,coilpath); % call the MARIE_COILLoad GUI
uiwait(gcf);

global COIL; % get the RHBM as global

% call the function to update the Geometry Plot
str = get(handles.Plot_ONOFF_Button,'String');    % get the current pushbutton string 
Plot_on = strcmp(str,'Plot ON'); 
Azimut = get(handles.Azimut_Slider,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
update_Geometry(handles.Geometry,Azimut,Elevation,Plot_on);

% update the entries of the model
if ~isempty(COIL.name)
    set(handles.ModelName_COIL,'String',COIL.name);
end
% Commit the new struct element to appdata
guidata(hObject, handles); 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Create_COIL.
function Create_COIL_Callback(hObject, eventdata, handles)
% hObject    handle to Create_COIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ispc
   % define the path to the coil
    coilpath = strcat(handles.path, 'coils\');
else
    % define the path to the coil
    coilpath = strcat(handles.path, 'coils/');
end

MARIE_SCOILGenerate(handles.Ncoils,coilpath); % call the MARIE_SCOILgenerate GUI
uiwait(gcf);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function ModelName_COIL_Callback(hObject, eventdata, handles)
% hObject    handle to ModelName_COIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ModelName_COIL as text
%        str2double(get(hObject,'String')) returns contents of ModelName_COIL as a double


% --- Executes during object creation, after setting all properties.
function ModelName_COIL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModelName_COIL (see GCBO)
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
function ModelName_RHBM_Callback(hObject, eventdata, handles)
% hObject    handle to ModelName_RHBM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ModelName_RHBM as text
%        str2double(get(hObject,'String')) returns contents of ModelName_RHBM as a double


% --- Executes during object creation, after setting all properties.
function ModelName_RHBM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ModelName_RHBM (see GCBO)
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
% --- Executes on button press in Create_RHBM.
function Create_RHBM_Callback(hObject, eventdata, handles)
% hObject    handle to Create_RHBM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_RHBMGenerate; % call the MARIE_RHBMGenerate GUI
uiwait(gcf);

global RHBM; % get the RHBM as global

% call the function to update the Geometry Plot
str = get(handles.Plot_ONOFF_Button,'String');    % get the current pushbutton string 
Plot_on = strcmp(str,'Plot ON'); 
Azimut = get(handles.Azimut_Slider,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
update_Geometry(handles.Geometry,Azimut,Elevation,Plot_on);

% update the entries of the model
if ~isempty(RHBM.name)
    set(handles.ModelName_RHBM,'String',RHBM.name);
end
% Commit the new struct element to appdata
guidata(hObject, handles); 


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Select_RHBM.
function Select_RHBM_Callback(hObject, eventdata, handles)
% hObject    handle to Select_RHBM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ispc
   % define the path to the coil
    rhbmpath = strcat(handles.path, 'bodies\');
else
    % define the path to the coil
    rhbmpath = strcat(handles.path, 'bodies/');
end

MARIE_RHBMLoad(handles.Nbodies,rhbmpath); % call the MARIE_RHBMLoad GUI
uiwait(gcf);

global RHBM; % get the RHBM as global

% call the function to update the Geometry Plot
str = get(handles.Plot_ONOFF_Button,'String');    % get the current pushbutton string 
Plot_on = strcmp(str,'Plot ON'); 
Azimut = get(handles.Azimut_Slider,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
update_Geometry(handles.Geometry,Azimut,Elevation,Plot_on);

% update the entries of the model
if ~isempty(RHBM.name)
    set(handles.ModelName_RHBM,'String',RHBM.name);
end
% Commit the new struct element to appdata
guidata(hObject, handles); 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Modify_RHBM.
function Modify_RHBM_Callback(hObject, eventdata, handles)
% hObject    handle to Modify_RHBM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_RHBMModify(handles.Geometry); % call the modify GUI
uiwait(gcf);

% call the function to update the Geometry Plot
str = get(handles.Plot_ONOFF_Button,'String');    % get the current pushbutton string 
Plot_on = strcmp(str,'Plot ON'); 
Azimut = get(handles.Azimut_Slider,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
update_Geometry(handles.Geometry,Azimut,Elevation,Plot_on);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Properties_RHBM.
function Properties_RHBM_Callback(hObject, eventdata, handles)
% hObject    handle to Properties_RHBM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_RHBMProperties; % Show model properties

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Export_SnP.
function Export_SnP_Callback(hObject, eventdata, handles)
% hObject    handle to Export_SnP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SOL;

if isempty(SOL.Zparam)
    uiwait(msgbox('Error: There is no data available','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
    uiresume(gcbf);
else
    if (size(SOL.Zparam,3) ~= length(SOL.freq))
        uiwait(msgbox('Error: Data is no coherent','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
        uiresume(gcbf);
    else
        
        MARIE_ExportSnP; % call function
        uiwait(gcf);        
    end
end



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on slider movement.
function Elevation_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Elevation_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Elevation = get(hObject,'Value');
Azimut = get(handles.Azimut_Slider,'Value');
axes(handles.Geometry);
view(Azimut,Elevation);
 axis image;
drawnow;


% --- Executes during object creation, after setting all properties.
function Elevation_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Elevation_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on slider movement.
function Azimut_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Azimut_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

Azimut = get(hObject,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
axes(handles.Geometry);
view(Azimut,Elevation);
 axis image;
drawnow;


% --- Executes during object creation, after setting all properties.
function Azimut_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Azimut_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Plot_ONOFF_Button.
function Plot_ONOFF_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_ONOFF_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = get(hObject,'String');    % get the current pushbutton string 
% Determine which string the label matches 
state = find(strcmp(str,handles.Strings)); 
% Toggle the button label to other string 
set(hObject,'String',handles.Strings{3-state}); 


if (state == 1) % If the index when entering was 1, Plot is on, so turn it off
    
    % clear the figure
    Azimut = get(handles.Azimut_Slider,'Value');
    Elevation = get(handles.Elevation_Slider,'Value');
    update_Geometry(handles.Geometry,Azimut,Elevation,0);
    
else % Plot is off, so turn it on
    
    % plot the figure
    Azimut = get(handles.Azimut_Slider,'Value');
    Elevation = get(handles.Elevation_Slider,'Value');
    update_Geometry(handles.Geometry,Azimut,Elevation,1);
    
 end


% --- Executes during object creation, after setting all properties.
function Plot_ONOFF_Button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot_ONOFF_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

handles.Strings = {'Plot ON';'Plot OFF'}; 
% Commit the new struct element to appdata
guidata(hObject, handles); 

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% function to update the geometry plot when a model is loaded
function update_Geometry(Geoaxes,Azimut,Elevation,Plot_on)

% get global variables with data
global RHBM;
global COIL;

% define view
axes(Geoaxes);
cla(Geoaxes);
view(Azimut,Elevation);
hold on;

xmin = []; xmax = [];
ymin = []; ymax = [];
zmin = []; zmax = [];

if (Plot_on)
    % plot RHBM
    if ~isempty(RHBM.r)
        
        mask_duke = 0*RHBM.epsilon_r;
        idxS = find( abs(RHBM.epsilon_r(:) - 1 + RHBM.sigma_e(:)));
        mask_duke(idxS) = 1;
        x = RHBM.r(:,:,:,1); xmin = min(x(:)); xmax = max(x(:));
        y = RHBM.r(:,:,:,2); ymin = min(y(:)); ymax = max(y(:));
        z = RHBM.r(:,:,:,3); zmin = min(z(:)); zmax = max(z(:));
        FV = isosurface(x,y,z,mask_duke);
        patch(FV, 'FaceColor','b', 'EdgeColor', 'none');
        hold on;
        clear x; clear y; clear z;
        
    end
    
    % plot COIL model
    if ~isempty(COIL.type)
        
        if(COIL.type == 'S')
            
            for ii = 1:size(COIL.elem,2)
                
                clear x; clear y; clear z;
                p1 = COIL.elem(1,ii);
                p2 = COIL.elem(2,ii);
                p3 = COIL.elem(3,ii);
                x(1) = COIL.node(1,p1); x(2) = COIL.node(1,p2); x(3) = COIL.node(1,p3);
                y(1) = COIL.node(2,p1); y(2) = COIL.node(2,p2); y(3) = COIL.node(2,p3);
                z(1) = COIL.node(3,p1); z(2) = COIL.node(3,p2); z(3) = COIL.node(3,p3);

                hold on;
                orange = [1 0.5 0.2];
                patch(x,y,z,orange, 'EdgeColor', 'none');
                xmin = min([xmin; x(:)]); xmax = max([xmax; x(:)]);
                ymin = min([ymin; y(:)]); ymax = max([ymax; y(:)]);
                zmin = min([zmin; z(:)]); zmax = max([zmax; z(:)]);
                clear x; clear y; clear z;
            end
            
        end
        
        if(COIL.type == 'W')
            
            for ii = 1:size(COIL.Pcoil,1)
                
                clear x; clear y; clear z;
                x(1) = COIL.Ncoil(ii,1); x(2) = COIL.Pcoil(ii,1);
                y(1) = COIL.Ncoil(ii,2); y(2) = COIL.Pcoil(ii,2);
                z(1) = COIL.Ncoil(ii,3); z(2) = COIL.Pcoil(ii,3);
                
                hold on;
                orange = [1 0.5 0.2];
                plot3(x,y,z,'Color', orange, 'LineWidth', 3.0);
                xmin = min([xmin; x(:)]); xmax = max([xmax; x(:)]);
                ymin = min([ymin; y(:)]); ymax = max([ymax; y(:)]);
                zmin = min([zmin; z(:)]); zmax = max([zmax; z(:)]);
                clear x; clear y; clear z;
                
            end
            
        end
        
    end
    
end

% set plot limits and stuff
if isempty(xmin)
    xmin = -1;
end
if isempty(xmax)
    xmax = 1;
else
    if xmax == xmin
        xmax = xmin + 0.001;
    end
end
if isempty(ymin)
    ymin = -1;
end
if isempty(ymax)
    ymax = 1;
else
    if ymax == ymin
        ymax = ymin + 0.001;
    end
end
if isempty(zmin)
    zmin = -1;
end
if isempty(zmax)
    zmax = 1;
else
    if zmax == zmin
        zmax = zmin + 0.001;
    end
end

set(gca, 'XLim', [xmin xmax]); xlabel(gca,'x (m)','Fontsize', 14);
set(gca, 'YLim', [ymin ymax]); ylabel(gca,'y (m)','Fontsize', 14);
set(gca, 'ZLim', [zmin zmax]); zlabel(gca,'z (m)','Fontsize', 14);
set(gca, 'DataAspectRatio', [1,1,1]);
grid on;
camlight;
lighting phong;
light('Position',[1 1 1],'Style','infinite');
light('Position',[-1 -1 1],'Style','infinite');
hold off;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RHBM;
global COIL;
global SOL;
global PULSE;

% get file name
[status_filename, status_path] = uiputfile('*.mat', 'Save current status as');

filename = sprintf('%s%s', status_path, status_filename);

% check that return was correct
if (status_filename(1) ~= 0) && (status_path(1) ~= 0)
    save(filename, 'RHBM', 'COIL', 'SOL', 'PULSE', '-v7.3');
    uiwait(msgbox('Current Data Saved','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
    uiresume(gcbf);
else
    f = warndlg('Incorrect file selected', 'Wrong File');
    waitfor(f);
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% get file name
[status_filename, status_path] = uigetfile('*.mat', 'Load status from');

filename = sprintf('%s%s', status_path, status_filename);
% check that return was correct
if (status_filename(1) ~= 0) && (status_path(1) ~= 0)
    
    clear RHBM;
    clear COIL;
    clear SOL;
    clear PULSE;
    
    global RHBM;
    global COIL;
    global SOL;
    global PULSE;
    
    load(filename);
    
%     if ~exist(RHBM.name)               
%         RHBM.name = 'No Model Selected';
%         RHBM.r = [];
%         RHBM.epsilon_r = [];
%         RHBM.sigma_e = [];
%         RHBM.rho = [];
%         RHBM.idxS = [];
%     end
%     
%     if ~exist(COIL.name)
%         COIL.name = 'No Model Selected';
%         COIL.type = 'N';
%     end
%     
%     if ~exist(SOL.Jsol)
%         SOL.Jsol = [];
%         SOL.Esol = [];
%         SOL.Bsol = [];
%         SOL.Ssol = [];
%         SOL.Gsar = [];
%         SOL.Pabs = [];
%         SOL.Sparam = [];
%         SOL.freq = [];
%     end

    uiwait(msgbox('Data Loaded','MARIE', 'custom', imread('Mlogo_blueorange','jpeg')));
    uiresume(gcbf);

else
    f = warndlg('Incorrect file selected', 'Wrong File');
    waitfor(f);
end
    
% call the function to update the Geometry Plot
str = get(handles.Plot_ONOFF_Button,'String');    % get the current pushbutton string 
Plot_on = strcmp(str,'Plot ON'); 
Azimut = get(handles.Azimut_Slider,'Value');
Elevation = get(handles.Elevation_Slider,'Value');
update_Geometry(handles.Geometry,Azimut,Elevation,Plot_on);

% update the entries of the models
if ~isempty(RHBM.name)
    set(handles.ModelName_RHBM,'String',RHBM.name);
end
if ~isempty(COIL.name)
    set(handles.ModelName_COIL,'String',COIL.name);
end
% Commit the new struct element to appdata
guidata(hObject, handles); 


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Apply_Pulse.
function Apply_Pulse_Callback(hObject, eventdata, handles)
% hObject    handle to Apply_Pulse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MARIE_Pulse;
uiwait(gcf);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
