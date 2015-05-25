function varargout = MARIE_Visualize(varargin)
% MARIE_VISUALIZE MATLAB code for MARIE_Visualize.fig
%      MARIE_VISUALIZE, by itself, creates a new MARIE_VISUALIZE or raises the existing
%      singleton*.
%
%      H = MARIE_VISUALIZE returns the handle to a new MARIE_VISUALIZE or the handle to
%      the existing singleton*.
%
%      MARIE_VISUALIZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_VISUALIZE.M with the given input arguments.
%
%      MARIE_VISUALIZE('Property','Value',...) creates a new MARIE_VISUALIZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_Visualize_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_Visualize_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_Visualize

% Last Modified by GUIDE v2.5 25-Jun-2014 19:37:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_Visualize_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_Visualize_OutputFcn, ...
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


% --- Executes just before MARIE_Visualize is made visible.
function MARIE_Visualize_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_Visualize (see VARARGIN)

% Choose default command line output for MARIE_Visualize
handles.output = hObject;

% initialize handles fields
handles.plottype = 'Plots';
handles.maptype = 'Maps';
handles.fieldresult = 'Available Results';
handles.freqresult = 'Available Results';
handles.mapview = 'Axial';


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_Visualize wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_Visualize_OutputFcn(hObject, eventdata, handles) 
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
% --- Executes on selection change in Freq_Results_Popup.
function Freq_Results_Popup_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_Results_Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Freq_Results_Popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Freq_Results_Popup

result = get(hObject,'Value');
freqresults = cellstr(get(hObject,'String'));
handles.freqresult = freqresults{result};
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Freq_Results_Popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq_Results_Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% generate the possible maps available as resuts
global SOL;
Nsol = 1;
freqresults{Nsol} = 'Available Results';

if ~isempty(SOL.Zparam)
    Nsol = Nsol+1;
    freqresults{Nsol} = 'S parameters';
    Nsol = Nsol+1;
    freqresults{Nsol} = 'Z parameters';
    Nsol = Nsol+1;
    freqresults{Nsol} = 'Y parameters';
end
if ~isempty(SOL.Pabs)
    Nsol = Nsol+1;
    freqresults{Nsol} = 'Absorbed Power';
end
if ~isempty(SOL.Gsar)
    Nsol = Nsol+1;
    freqresults{Nsol} = 'Global SAR';
end

set(hObject,'String',freqresults);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Freq_Plot_Button.
function Freq_Plot_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_Plot_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global SOL;
resulttype = get(handles.Freq_Results_Popup, 'Value');

switch resulttype
    
    case 1
        type = 'W';
    case 5
        vec = SOL.Pabs;
        freq = SOL.freq;
        type = 'P';
    case 6
        vec = SOL.Gsar;
        freq = SOL.freq;
        type = 'G';
    otherwise
        vec = SOL.Zparam;
        freq = SOL.freq;
        type = 'S';
end

if (type == 'W')
    msgbox('Analysis not available for the selected result','Unavailable analysis','warn');
else
    Nports = size(vec,1);
    channelidx = 1:Nports;
    MARIE_PlotFreq(vec,freq,type,channelidx);
    uiwait(gcf);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on selection change in Field_Results_Popup.
function Field_Results_Popup_Callback(hObject, eventdata, handles)
% hObject    handle to Field_Results_Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Field_Results_Popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Field_Results_Popup

result = get(hObject,'Value');
fieldresults = cellstr(get(hObject,'String'));
handles.fieldresult = fieldresults{result}; 
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Field_Results_Popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Field_Results_Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% generate the possible maps available as resuts
global SOL;
global RHBM;
Nsol = 1;
fieldresults{Nsol} = 'Available Results';

if ~isempty(RHBM.epsilon_r)
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'Relative Permittivity (epsilon_r)';
end
if ~isempty(RHBM.sigma_e)
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'Electric Conductivity (sigma_e)';
end
if ~isempty(RHBM.rho)
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'Density (rho)';
end
if ~isempty(SOL.Jsol)
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'Electric Current (J)';
end
if ~isempty(SOL.Esol)
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'Electric Field (E)';
end
if ~isempty(SOL.Bsol)
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'Magnetic Field (B)';
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'B1+ component';
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'B1- component';
end
if ~isempty(SOL.Ssol)
    Nsol = Nsol+1;
    fieldresults{Nsol} = 'Local SAR (unweighted)';
end

set(hObject,'String',fieldresults);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on selection change in Map_Type_Popup.
function Map_Type_Popup_Callback(hObject, eventdata, handles)
% hObject    handle to Map_Type_Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Map_Type_Popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Map_Type_Popup

type = get(hObject,'Value');
maptypes = cellstr(get(hObject,'String'));
handles.maptype = maptypes{type};
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Map_Type_Popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Map_Type_Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% generate the possible maps available as resuts
maptypes = {'Maps', 'Magnitude', 'Phase', 'RMS', 'X comp. Magnitude', 'X comp. Phase', 'Y comp. Magnitude', 'Y comp. Phase', 'Z comp. Magnitude', 'Z comp. Phase'};
set(hObject,'String',maptypes);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on selection change in Map_View_Popup.
function Map_View_Popup_Callback(hObject, eventdata, handles)
% hObject    handle to Map_View_Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Map_View_Popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Map_View_Popup

view = get(hObject,'Value');
mapviews = cellstr(get(hObject,'String'));
handles.mapview = mapviews{view};
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Map_View_Popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Map_View_Popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% generate the possible maps available as resuts
mapviews = {'Axial', 'Saggital', 'Coronal'};
set(hObject,'String',mapviews);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Maps_Plot_Button.
function Maps_Plot_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Maps_Plot_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% call the plot for the desired map
global RHBM;
global SOL;


% select the mapview
mapview = 'A';
if (strcmp(handles.mapview, 'Saggital'))
    mapview = 'S';
end
if (strcmp(handles.mapview, 'Coronal'))
    mapview = 'C';
end

% get the number related to the map type
maptypenum = find(strcmp(handles.maptype, cellstr(get(handles.Map_Type_Popup,'String'))));


% select the result to plot
vec = []; wrongtype = 0; Nports = 0; Nfreqs = 0;
if strcmp(handles.fieldresult, 'Relative Permittivity (epsilon_r)')
    switch maptypenum
        case 2 %'Magnitude'
            nameplot = strcat(RHBM.name, ' -- epsilon_r');
            vec = RHBM.epsilon_r;
        otherwise
            wrongtype = 1;
    end
end
if strcmp(handles.fieldresult, 'Electric Conductivity (sigma_e)')
    switch maptypenum
        case 2 %'Magnitude'
            nameplot = strcat(RHBM.name, ' -- sigma_e');
            vec = RHBM.sigma_e;
        otherwise
            wrongtype = 1;
    end
end
if strcmp(handles.fieldresult, 'Density (rho)')
    switch maptypenum
        case 2 %'Magnitude'
            nameplot = strcat(RHBM.name, ' -- rho');
            vec = RHBM.rho;
        otherwise
            wrongtype = 1;
    end
end


if strcmp(handles.fieldresult, 'Electric Current (J)')
    [L,M,N,~,Nports,Nfreqs] = size(SOL.Jsol);
    solvec = SOL.Jsol;
    vec = zeros(L,M,N,Nports,Nfreqs);
    switch maptypenum
        case 2 %'Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(J)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(sum(solvec(:,:,:,:,ii,jj),4));
                end
            end
        case 3 %'Phase'
            nameplot = strcat(RHBM.name, ' -- phase(J)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(sum(solvec(:,:,:,:,ii,jj),4));
                end
            end
        case 4 %'RMS'
            nameplot = strcat(RHBM.name, ' -- RMS(J)');
            solvec = solvec.*conj(solvec);
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = sqrt(sum(solvec(:,:,:,:,ii,jj),4)./2);
                end
            end
        case 5 % 'X comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(J_x)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        case 6 %'X comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(J_x)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        case 7 % 'Y comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(J_y)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,2,ii,jj)));
                end
            end
        case 8 %'Y comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(J_y)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,2,ii,jj)));
                end
            end         
        case 9 % 'Z comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(J_z)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,3,ii,jj)));
                end
            end
        case 10 %'Z comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(J_z)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,3,ii,jj)));
                end
            end
        otherwise
            wrongtype = 1;
    end
end

if strcmp(handles.fieldresult, 'Electric Field (E)')
    [L,M,N,~,Nports,Nfreqs] = size(SOL.Esol);
    solvec = SOL.Esol;
    vec = zeros(L,M,N,Nports,Nfreqs);
    switch maptypenum
        case 2 %'Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(E)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(sum(solvec(:,:,:,:,ii,jj),4));
                end
            end
        case 3 %'Phase'
            nameplot = strcat(RHBM.name, ' -- phase(E)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(sum(solvec(:,:,:,:,ii,jj),4));
                end
            end
        case 4 %'RMS'
            nameplot = strcat(RHBM.name, ' -- RMS(E)');
            solvec = solvec.*conj(solvec);
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = sqrt(sum(solvec(:,:,:,:,ii,jj),4)./2);
                end
            end
        case 5 % 'X comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(E_x)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        case 6 %'X comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(E_x)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        case 7 % 'Y comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(E_y)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,2,ii,jj)));
                end
            end
        case 8 %'Y comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(E_y)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,2,ii,jj)));
                end
            end         
        case 9 % 'Z comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(E_z)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,3,ii,jj)));
                end
            end
        case 10 %'Z comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(E_z)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,3,ii,jj)));
                end
            end
        otherwise
            wrongtype = 1;
    end
end

if strcmp(handles.fieldresult, 'Magnetic Field (B)')
    [L,M,N,~,Nports,Nfreqs] = size(SOL.Bsol);
    solvec = SOL.Bsol;
    vec = zeros(L,M,N,Nports,Nfreqs);
    switch maptypenum
        case 2 %'Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(B)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(sum(solvec(:,:,:,:,ii,jj),4));
                end
            end
        case 3 %'Phase'
            nameplot = strcat(RHBM.name, ' -- phase(B)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(sum(solvec(:,:,:,:,ii,jj),4));
                end
            end
        case 4 %'RMS'
            nameplot = strcat(RHBM.name, ' -- RMS(B)');
            solvec = solvec.*conj(solvec);
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = sqrt(sum(solvec(:,:,:,:,ii,jj),4)./2);
                end
            end
        case 5 % 'X comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(B_x)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        case 6 %'X comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(B_x)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        case 7 % 'Y comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(B_y)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,2,ii,jj)));
                end
            end
        case 8 %'Y comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(B_y)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,2,ii,jj)));
                end
            end         
        case 9 % 'Z comp. Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(B_z)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,3,ii,jj)));
                end
            end
        case 10 %'Z comp. Phase'
            nameplot = strcat(RHBM.name, ' -- phase(B_z)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,3,ii,jj)));
                end
            end
        otherwise
            wrongtype = 1;
    end
end
if strcmp(handles.fieldresult, 'B1+ component')
    [L,M,N,~,Nports,Nfreqs] = size(SOL.Bsol);
    vec = zeros(L,M,N,Nports,Nfreqs);
    solvec = (SOL.Bsol(:,:,:,1,:,:) + 1j*SOL.Bsol(:,:,:,2,:,:))/2;
    switch maptypenum
        case 2 %'Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(B1+)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        case 3 %'Phase'
            nameplot = strcat(RHBM.name, ' -- phase(B1+)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        otherwise
            wrongtype = 1;
    end
end
if strcmp(handles.fieldresult, 'B1- component')
    [L,M,N,~,Nports,Nfreqs] = size(SOL.Bsol);
    vec = zeros(L,M,N,Nports,Nfreqs);
    solvec = (SOL.Bsol(:,:,:,1) - 1j*SOL.Bsol(:,:,:,2))/2;
    switch maptypenum
        case 2 %'Magnitude'
            nameplot = strcat(RHBM.name, ' -- abs(B1-)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = abs(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        case 3 %'Phase'
            nameplot = strcat(RHBM.name, ' -- phase(B1-)');
            for ii = 1:Nports
                for jj = 1:Nfreqs
                    vec(:,:,:,ii,jj) = angle(squeeze(solvec(:,:,:,1,ii,jj)));
                end
            end
        otherwise
            wrongtype = 1;
    end
end
if strcmp(handles.fieldresult, 'Local SAR (unweighted)')
    switch maptypenum
        case 2 %'Magnitude'
            [~,~,~,Nports,Nfreqs] = size(SOL.Ssol);
            nameplot = strcat(RHBM.name, ' -- Local SAR');
            vec = SOL.Ssol;
        otherwise
            wrongtype = 1;
    end
end

if wrongtype
    msgbox('Analysis not available for the selected result','Unavailable analysis','warn');
else
    if isempty(vec)
        msgbox('Empty result vector: select vector','Unavailable vector','warn');
    else
        if (Nports > 1) || (Nfreqs > 1)
            % call the plot
            channelidx = 1:Nports;
            MARIE_PlotSlice_Multi(vec,RHBM.r,mapview,nameplot,RHBM.idxS,SOL.freq,channelidx);
            uiwait(gcf);
        else
            vec = squeeze(vec);
            MARIE_PlotSlice(vec,RHBM.r,mapview,nameplot,RHBM.idxS);
            uiwait(gcf);
        end
    end
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
