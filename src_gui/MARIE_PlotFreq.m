function varargout = MARIE_PlotFreq(varargin)
% MARIE_PLOTFREQ MATLAB code for MARIE_PlotFreq.fig
%      MARIE_PLOTFREQ, by itself, creates a new MARIE_PLOTFREQ or raises the existing
%      singleton*.
%
%      H = MARIE_PLOTFREQ returns the handle to a new MARIE_PLOTFREQ or the handle to
%      the existing singleton*.
%
%      MARIE_PLOTFREQ('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_PLOTFREQ.M with the given input arguments.
%
%      MARIE_PLOTFREQ('Property','Value',...) creates a new MARIE_PLOTFREQ or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_PlotFreq_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_PlotFreq_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_PlotFreq

% Last Modified by GUIDE v2.5 25-Jun-2014 19:33:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_PlotFreq_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_PlotFreq_OutputFcn, ...
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


% --- Executes just before MARIE_PlotFreq is made visible.
function MARIE_PlotFreq_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_PlotFreq (see VARARGIN)

% Choose default command line output for MARIE_PlotFreq
handles.output = hObject;

% get the variables
handles.vector = varargin{1};
handles.freqval = varargin{2};
handles.type = varargin{3};
handles.channel = varargin{4};

handles.plotlist = []; % curves to plot
handles.namelist = []; % names of the curves

switch handles.type
 case 'S' % s-parameter plots
     
     % set channel values
     Nin = size(handles.vector,1);
     for ii = 1:Nin
         Ninvals{ii} = sprintf('%g', handles.channel(ii));
     end
     set(handles.In_Box_Select, 'String', Ninvals);
     set(handles.In_Box_Select, 'Value', 1);
     set(handles.Out_Box_Select, 'String', Ninvals);
     set(handles.Out_Box_Select, 'Value', 1);

     plottypes = {'Magnitude', 'Phase', 'Real', 'Imaginary', 'Magnitude (dB)'};
     set(handles.Plot_Add_Box_Select, 'String', plottypes);
     set(handles.Plot_Add_Box_Select, 'Value', 1);

     resulttypes = {'Z parameter', 'Y parameter', 'S parameter'};
     set(handles.Result_Add_Box_Select, 'String', resulttypes);
     set(handles.Result_Add_Box_Select, 'Value', 1);

 case 'P'

     % set channel values
     Nin = size(handles.vector,1);
     for ii = 1:Nin
         Ninvals{ii} = sprintf('%g', handles.channel(ii));
     end
     set(handles.In_Box_Select, 'String', Ninvals);
     set(handles.In_Box_Select, 'Value', 1);

     plottypes = {'Magnitude', 'Magnitude (dB)'};
     set(handles.Plot_Add_Box_Select, 'String', plottypes);
     set(handles.Plot_Add_Box_Select, 'Value', 1);

     resulttypes = {'Power absorbed'};
     set(handles.Result_Add_Box_Select, 'String', resulttypes);
     set(handles.Result_Add_Box_Select, 'Value', 1);

 otherwise

     % set channel values
     Nin = size(handles.vector,1);
     for ii = 1:Nin
         Ninvals{ii} = sprintf('%g', handles.channel(ii));
     end
     set(handles.In_Box_Select, 'String', Ninvals);
     set(handles.In_Box_Select, 'Value', 1);

     plottypes = {'Magnitude', 'Magnitude (dB)'};
     set(handles.Plot_Add_Box_Select, 'String', plottypes);
     set(handles.Plot_Add_Box_Select, 'Value', 1);

     resulttypes = {'Global SAR'};
     set(handles.Result_Add_Box_Select, 'String', resulttypes);
     set(handles.Result_Add_Box_Select, 'Value', 1);

end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_PlotFreq wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% --- Outputs from this function are returned to the command line.
function varargout = MARIE_PlotFreq_OutputFcn(hObject, eventdata, handles) 
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
% --- Executes on button press in Remove_Button.
function Remove_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

remove_number = get(handles.Result_Remove_Box_Select, 'Value');
plot_list = cellstr(get(handles.Result_Remove_Box_Select, 'String'));

% remove plot from list
newplotlist = [];
newnamelist = [];
Nplots = length(handles.plotlist);
newNplots = 0;
for ii = 1: Nplots
    if ~strcmp(handles.namelist{ii},plot_list{remove_number})
        newNplots = newNplots+1;
        newplotlist{newNplots} = handles.plotlist{ii};
        newnamelist{newNplots} = handles.namelist{ii};
    end
end

% add the plot to the current list
handles.plotlist = newplotlist;
handles.namelist = newnamelist;

% add the plot to the current list: able to be removed
set(handles.Result_Remove_Box_Select, 'String', newnamelist);
set(handles.Result_Remove_Box_Select, 'Value', 1);

guidata(hObject, handles);

% generate new plot
GeneratePlot(handles);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% --- Executes on selection change in Result_Remove_Box_Select.
function Result_Remove_Box_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Result_Remove_Box_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Result_Remove_Box_Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Result_Remove_Box_Select


% --- Executes during object creation, after setting all properties.
function Result_Remove_Box_Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Result_Remove_Box_Select (see GCBO)
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
% --- Executes on selection change in In_Box_Select.
function In_Box_Select_Callback(hObject, eventdata, handles)
% hObject    handle to In_Box_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns In_Box_Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from In_Box_Select


% --- Executes during object creation, after setting all properties.
function In_Box_Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to In_Box_Select (see GCBO)
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
% --- Executes on selection change in Out_Box_Select.
function Out_Box_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Out_Box_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Out_Box_Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Out_Box_Select


% --- Executes during object creation, after setting all properties.
function Out_Box_Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Out_Box_Select (see GCBO)
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
% --- Executes on selection change in Result_Add_Box_Select.
function Result_Add_Box_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Result_Add_Box_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Result_Add_Box_Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Result_Add_Box_Select


% --- Executes during object creation, after setting all properties.
function Result_Add_Box_Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Result_Add_Box_Select (see GCBO)
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
% --- Executes on selection change in Plot_Add_Box_Select.
function Plot_Add_Box_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_Add_Box_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Plot_Add_Box_Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Plot_Add_Box_Select


% --- Executes during object creation, after setting all properties.
function Plot_Add_Box_Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Plot_Add_Box_Select (see GCBO)
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
% --- Executes on button press in Add_Button.
function Add_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Add_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the data to plot
resultnumber = get(handles.Result_Add_Box_Select, 'Value');
plotnumber = get(handles.Plot_Add_Box_Select, 'Value');
plotname = cellstr(get(handles.Plot_Add_Box_Select, 'String'));
Nin = get(handles.In_Box_Select, 'Value');
Nout = get(handles.Out_Box_Select, 'Value');

Nplots = length(handles.plotlist);

% get the type of plot
switch handles.type
 case 'S' % s-parameter plots
     
     switch resultnumber
         case 1 % Z-parameters
             vec = handles.vector;
             name = 'Z';
         case 2 % Y-parameters
             vec = z2y(handles.vector);
             name = 'Y';
         case 3 % S parameters
             vec = z2s(handles.vector,50);
             name = 'S';
     end

     switch plotnumber
         case 1 % magnitude
             vec = squeeze(abs(vec(Nout,Nin,:)));
         case 2 % phase
             vec = squeeze(angle(vec(Nout,Nin,:)));
         case 3 % real
             vec = squeeze(real(vec(Nout,Nin,:)));
         case 4 % imag
             vec = squeeze(imag(vec(Nout,Nin,:)));
         case 5 % dB
             vec = squeeze(20*log10(abs(vec(Nout,Nin,:))));
     end

     handles.namelist{Nplots+1} = sprintf('%s_%s_%d_%d',plotname{plotnumber},name,Nout,Nin);
     handles.plotlist{Nplots+1} = vec;

 case 'P'

     vec = handles.vector;
     switch plotnumber
         case 1 % magnitude
             vec = squeeze(abs(vec(Nin,:)));
         case 2 % dB
             vec = squeeze(20*log10(abs(vec(Nin,:))));
     end

     handles.namelist{Nplots+1} = sprintf('%s Pabs_ch%d',plotname{plotnumber},Nin);
     handles.plotlist{Nplots+1} = vec;

 otherwise

     vec = handles.vector;
     switch plotnumber
         case 1 % magnitude
             vec = squeeze(abs(vec(Nin,:)));
         case 2 % dB
             vec = squeeze(20*log10(abs(vec(Nin,:))));
     end

     handles.namelist{Nplots+1} = sprintf('%s GSAR_ch%d',plotname{plotnumber},Nin);
     handles.plotlist{Nplots+1} = vec;

end

% add the plot to the current list: able to be removed
set(handles.Result_Remove_Box_Select, 'String', handles.namelist);
set(handles.Result_Remove_Box_Select, 'Value', 1);

% update handles
guidata(hObject, handles);

% generate new plot
GeneratePlot(handles);


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% function to plot
function GeneratePlot(handles)

% get data
currentaxes = handles.axes1;
% clear current axes
cla(currentaxes);

ff = squeeze(handles.freqval/1e6);

% plot
Nplots = length(handles.plotlist);
for ii = 1:Nplots
        
    vec = squeeze(handles.plotlist{ii});
    plot(ff,vec,'-', 'Color', rand(3,1), 'LineWidth', 2);
    hold on;
        
end

hold off;
legend(handles.namelist);
ylabel(gca,'Value','Fontsize', 14);
xlabel(gca,'Frequency (MHz)','Fontsize', 14);
set(gca, 'FontSize', 18);
grid on;


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
figure(FIGIDX);


ff = squeeze(handles.freqval/1e6);

% plot
Nplots = length(handles.plotlist);
for ii = 1:Nplots
        
    vec = squeeze(handles.plotlist{ii});
    plot(ff,vec,'-', 'Color', rand(3,1), 'LineWidth', 2);
    hold on;
        
end

hold off;
legend(handles.namelist);
ylabel(gca,'Value','Fontsize', 14);
xlabel(gca,'Frequency (MHz)','Fontsize', 14);
set(gca, 'FontSize', 18);
grid on;

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
