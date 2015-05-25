function varargout = MARIE_PlotSlice_Multi(varargin)
% MARIE_PLOTSLICE_MULTI MATLAB code for MARIE_PlotSlice_Multi.fig
%      MARIE_PLOTSLICE_MULTI, by itself, creates a new MARIE_PLOTSLICE_MULTI or raises the existing
%      singleton*.
%
%      H = MARIE_PLOTSLICE_MULTI returns the handle to a new MARIE_PLOTSLICE_MULTI or the handle to
%      the existing singleton*.
%
%      MARIE_PLOTSLICE_MULTI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_PLOTSLICE_MULTI.M with the given input arguments.
%
%      MARIE_PLOTSLICE_MULTI('Property','Value',...) creates a new MARIE_PLOTSLICE_MULTI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_PlotSlice_Multi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_PlotSlice_Multi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_PlotSlice_Multi

% Last Modified by GUIDE v2.5 25-Jun-2014 19:44:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_PlotSlice_Multi_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_PlotSlice_Multi_OutputFcn, ...
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


% --- Executes just before MARIE_PlotSlice_Multi is made visible.
function MARIE_PlotSlice_Multi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_PlotSlice_Multi (see VARARGIN)

% Choose default command line output for MARIE_PlotSlice_Multi
handles.output = hObject;

% get the variables
handles.vector = varargin{1}; % vector to plot
handles.r = varargin{2}; % domain
handles.sliceview = varargin{3}; % view of the cut
handles.figurename = varargin{4}; % name of the map
handles.idxB = varargin{5}; % positions of the body
handles.freqs = varargin{6}; % frequencies vector
handles.channel = varargin{7}; % channels or excitation
[L,M,N,~] = size(handles.r);
nD = L*M*N;
handles.idxD = squeeze(1:nD);
handles.idxA = setdiff(handles.idxD, handles.idxB);

% find values of domain
xvec = unique(handles.r(:,1,1,1));
yvec = unique(handles.r(1,:,1,2));
zvec = unique(handles.r(1,1,:,3));

% set limits for the slide bars and define title
sliceview = handles.sliceview;
switch sliceview
    case 'S'        
        % -------------------------------------------
        % Sag
        figurename = '(Sag.) -- ';
        % set limits to slice bar
        set(handles.Pos_Slider,'Max',max(xvec(:)));
        set(handles.Pos_Slider,'Min',min(xvec(:)));
        slicepos = min(abs(xvec(:)));
        % set value of the Slide position
        set(handles.Pos_Box, 'String', sprintf('x = %g', slicepos));
        
    case 'C'
        % -------------------------------------------
        % Cor
        figurename = '(Cor.) -- ';
        % set limits to slice bar
        set(handles.Pos_Slider,'Max',max(yvec(:)));
        set(handles.Pos_Slider,'Min',min(yvec(:)));
        slicepos = min(abs(yvec(:)));
        % set value of the Slide position
        set(handles.Pos_Box, 'String', sprintf('y = %g', slicepos));
        
    otherwise % axial
        % -------------------------------------------
        % Axial
        figurename = '(Axial) -- ';
        % set limits to slice bar
        set(handles.Pos_Slider,'Max',max(zvec(:)));
        set(handles.Pos_Slider,'Min',min(zvec(:)));
        slicepos = min(abs(zvec(:)));
        % set value of the Slide position
        set(handles.Pos_Box, 'String', sprintf('z = %g', slicepos));
        
end

% set title
figurename = strcat(figurename, handles.figurename);
set(handles.Figure_Title,'String',figurename);

% set frequency values
for ii = 1:length(handles.freqs(:))
    freqvals{ii} = sprintf('%g Hz',handles.freqs(ii));
end
set(handles.Freq_Values_Menu, 'String', freqvals);
set(handles.Freq_Values_Menu, 'Value', 1);

% set channel values
Nchannels = size(handles.vector,4);
for ii = 1:Nchannels
    channelvals{ii} = sprintf('%g', handles.channel(ii));
end
set(handles.Channel_Values_Menu, 'String', channelvals);
set(handles.Channel_Values_Menu, 'Value', 1);

% call the plot generator
GenerateMap(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_PlotSlice_Multi wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% --- Outputs from this function are returned to the command line.
function varargout = MARIE_PlotSlice_Multi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on slider movement.
function Pos_Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% call the plot generator
GenerateMap(handles);


% --- Executes during object creation, after setting all properties.
function Pos_Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_Slider (see GCBO)
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
function Pos_Box_Callback(hObject, eventdata, handles)
% hObject    handle to Pos_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pos_Box as text
%        str2double(get(hObject,'String')) returns contents of Pos_Box as a double


% --- Executes during object creation, after setting all properties.
function Pos_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pos_Box (see GCBO)
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
function Figure_Title_Callback(hObject, eventdata, handles)
% hObject    handle to Figure_Title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Figure_Title as text
%        str2double(get(hObject,'String')) returns contents of Figure_Title as a double


% --- Executes during object creation, after setting all properties.
function Figure_Title_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Figure_Title (see GCBO)
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
% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on selection change in Freq_Values_Menu.
function Freq_Values_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Freq_Values_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Freq_Values_Menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Freq_Values_Menu

% call the plot generator
GenerateMap(handles);

% --- Executes during object creation, after setting all properties.
function Freq_Values_Menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq_Values_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on selection change in Channel_Values_Menu.
function Channel_Values_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Channel_Values_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Channel_Values_Menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Channel_Values_Menu

% call the plot generator
GenerateMap(handles);

% --- Executes during object creation, after setting all properties.
function Channel_Values_Menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Channel_Values_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Mask_Air_Select.
function Mask_Air_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Mask_Air_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mask_Air_Select

set(handles.Mask_All_Select,'Value',0);
set(handles.Mask_Body_Select,'Value',0);
% call plot generator
GenerateMap(handles);


% --- Executes on button press in Mask_Body_Select.
function Mask_Body_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Mask_Body_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mask_Body_Select

set(handles.Mask_Air_Select,'Value',0);
set(handles.Mask_All_Select,'Value',0);
% call plot generator
GenerateMap(handles);


% --- Executes on button press in Mask_All_Select.
function Mask_All_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Mask_All_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Mask_All_Select

set(handles.Mask_Air_Select,'Value',0);
set(handles.Mask_Body_Select,'Value',0);
% call plot generator
GenerateMap(handles);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% function to plot map
function GenerateMap(handles)

% get data
currentaxes = handles.axes1;
r = handles.r;
vector = handles.vector;
sliceview = handles.sliceview;

% get frequency position
freqpos = get(handles.Freq_Values_Menu,'Value');

% get the channel
channel = get(handles.Channel_Values_Menu,'Value');

% get the slice position
slicepos = get(handles.Pos_Slider,'Value');

% get the kind of mask
maskvec(1) = get(handles.Mask_All_Select,'Value');
maskvec(2) = get(handles.Mask_Air_Select,'Value');
maskvec(3) = get(handles.Mask_Body_Select,'Value');
maskval = find(maskvec);
switch maskval
    case 2
        idxMap = handles.idxA; % air
    case 3
        idxMap = handles.idxB; % body
    otherwise
        idxMap = handles.idxD; % whole domain
end


% get domain data
[L,M,N,~] = size(r);

% get corresponding vector, mask and scaling
vector = squeeze(vector(:,:,:,channel,freqpos));
vec = zeros(L,M,N);
vec(idxMap) = vector(idxMap);
scalemax = max(vec(:));
scalemin = min(vec(:));

% find indexes for the cuts
xvec = r(:,1,1,1);
yvec = r(1,:,1,2);
zvec = r(1,1,:,3);

switch sliceview
    case 'S'        
        % -------------------------------------------
        % Sag

        % set position in box
        [~, xidx] = min(abs(xvec - slicepos));
        xpos = xvec(xidx);
        set(handles.Pos_Box,'String',sprintf('x = %g',xpos));
        
        % plot to fix limits of scale
        Scaleplot = scalemax*ones(M,N);
        Scaleplot(end,end) = scalemin;
       
        % perform the cut and squeeze
        vecplot = squeeze(vec(xidx,:,:));
        vecplot = rot90(vecplot(end:-1:1,:)); %rotate to vertical and reverse
        
        % clear current axes
        cla(currentaxes);
        % plot
        imagesc(rot90(Scaleplot));
        hold on;
        imagesc(vecplot);
        xlabel(gca,'y','Fontsize', 14);
        ylabel(gca,'z','Fontsize', 14);
        colormap('hot'); colorbar;
        axis image; set(gca, 'FontSize', 18);
        
        
    case 'C'
        % -------------------------------------------
        % Cor
        
        % set position in box
        [~, yidx] = min(abs(yvec - slicepos));
        ypos = yvec(yidx);
        set(handles.Pos_Box,'String',sprintf('y = %g',ypos));

        % plot to fix limits of scale
        Scaleplot = scalemax*ones(L,N);
        Scaleplot(end,end) = scalemin;
       
        % perform the cut and squeeze
        vecplot = squeeze(vec(:,yidx,:));
        vecplot = rot90(vecplot); %rotate to vertical
        
        % clear current axes
        cla(currentaxes);
        % plot
        imagesc(rot90(Scaleplot));
        hold on;
        imagesc(vecplot);
        xlabel(gca,'x','Fontsize', 14);
        ylabel(gca,'z','Fontsize', 14);
        colormap('hot'); colorbar;
        axis image; set(gca, 'FontSize', 18);
        
        
    otherwise % axial
        % -------------------------------------------
        % Axial

        % set position in box
        [~, zidx] = min(abs(zvec - slicepos));
        zpos = zvec(zidx);
        set(handles.Pos_Box,'String',sprintf('z = %g',zpos));

        % plot to fix limits of scale
        Scaleplot = scalemax*ones(L,M);
        Scaleplot(end,end) = scalemin;

        % perform the cut and squeeze
        vecplot = squeeze(vec(:,:,zidx));
        vecplot = rot90(vecplot); %rotate to vertical
        
        % clear current axes
        cla(currentaxes);
        % plot
        imagesc(rot90(Scaleplot));
        hold on;
        imagesc(vecplot)
        xlabel(gca,'x','Fontsize', 14);
        ylabel(gca,'y','Fontsize', 14);
        colormap('hot'); colorbar;
        axis image; set(gca, 'FontSize', 18);
        
        
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on button press in Export_Button.
function Export_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Export_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global FIGIDX;
FIGIDX = FIGIDX + 1;


% get data
r = handles.r;
vector = handles.vector;
sliceview = handles.sliceview;

% get frequency position
freqpos = get(handles.Freq_Values_Menu,'Value');

% get the channel
channel = get(handles.Channel_Values_Menu,'Value');

% get the slice position
slicepos = get(handles.Pos_Slider,'Value');

% get the kind of mask
maskvec(1) = get(handles.Mask_All_Select,'Value');
maskvec(2) = get(handles.Mask_Air_Select,'Value');
maskvec(3) = get(handles.Mask_Body_Select,'Value');
maskval = find(maskvec);
switch maskval
    case 2
        idxMap = handles.idxA; % air
    case 3
        idxMap = handles.idxB; % body
    otherwise
        idxMap = handles.idxD; % whole domain
end


% get domain data
[L,M,N,~] = size(r);

% get corresponding vector, mask and scaling
vector = squeeze(vector(:,:,:,channel,freqpos));
vec = zeros(L,M,N);
vec(idxMap) = vector(idxMap);
scalemax = max(vec(:));
scalemin = min(vec(:));

% find indexes for the cuts
xvec = r(:,1,1,1);
yvec = r(1,:,1,2);
zvec = r(1,1,:,3);

switch sliceview
    case 'S'        
        % -------------------------------------------
        % Sag

        % set position in box
        [~, xidx] = min(abs(xvec - slicepos));
        xpos = xvec(xidx);
        set(handles.Pos_Box,'String',sprintf('x = %g',xpos));
        
        % plot to fix limits of scale
        Scaleplot = scalemax*ones(M,N);
        Scaleplot(end,end) = scalemin;
       
        % perform the cut and squeeze
        vecplot = squeeze(vec(xidx,:,:));
        vecplot = rot90(vecplot(end:-1:1,:)); %rotate to vertical and reverse
        
        % clear current axes
        figure(FIGIDX);
        % plot
        imagesc(rot90(Scaleplot));
        hold on;
        imagesc(vecplot);
        xlabel(gca,'y','Fontsize', 14);
        ylabel(gca,'z','Fontsize', 14);
        colormap('hot'); colorbar;
        axis image; set(gca, 'FontSize', 18);
        
        
    case 'C'
        % -------------------------------------------
        % Cor
        
        % set position in box
        [~, yidx] = min(abs(yvec - slicepos));
        ypos = yvec(yidx);
        set(handles.Pos_Box,'String',sprintf('y = %g',ypos));

        % plot to fix limits of scale
        Scaleplot = scalemax*ones(L,N);
        Scaleplot(end,end) = scalemin;
       
        % perform the cut and squeeze
        vecplot = squeeze(vec(:,yidx,:));
        vecplot = rot90(vecplot); %rotate to vertical
        
        % clear current axes
        figure(FIGIDX);
        % plot
        imagesc(rot90(Scaleplot));
        hold on;
        imagesc(vecplot);
        xlabel(gca,'x','Fontsize', 14);
        ylabel(gca,'z','Fontsize', 14);
        colormap('hot'); colorbar;
        axis image; set(gca, 'FontSize', 18);
        
        
    otherwise % axial
        % -------------------------------------------
        % Axial

        % set position in box
        [~, zidx] = min(abs(zvec - slicepos));
        zpos = zvec(zidx);
        set(handles.Pos_Box,'String',sprintf('z = %g',zpos));

        % plot to fix limits of scale
        Scaleplot = scalemax*ones(L,M);
        Scaleplot(end,end) = scalemin;

        % perform the cut and squeeze
        vecplot = squeeze(vec(:,:,zidx));
        vecplot = rot90(vecplot); %rotate to vertical
        
        % clear current axes
        figure(FIGIDX);
        % plot
        imagesc(rot90(Scaleplot));
        hold on;
        imagesc(vecplot)
        xlabel(gca,'x','Fontsize', 14);
        ylabel(gca,'y','Fontsize', 14);
        colormap('hot'); colorbar;
        axis image; set(gca, 'FontSize', 18);
        
        
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
