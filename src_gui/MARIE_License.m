function varargout = MARIE_License(varargin)
% MARIE_LICENSE MATLAB code for MARIE_License.fig
%      MARIE_LICENSE, by itself, creates a new MARIE_LICENSE or raises the existing
%      singleton*.
%
%      H = MARIE_LICENSE returns the handle to a new MARIE_LICENSE or the handle to
%      the existing singleton*.
%
%      MARIE_LICENSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MARIE_LICENSE.M with the given input arguments.
%
%      MARIE_LICENSE('Property','Value',...) creates a new MARIE_LICENSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MARIE_License_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MARIE_License_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MARIE_License

% Last Modified by GUIDE v2.5 24-Jun-2014 22:09:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MARIE_License_OpeningFcn, ...
                   'gui_OutputFcn',  @MARIE_License_OutputFcn, ...
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


% --- Executes just before MARIE_License is made visible.
function MARIE_License_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MARIE_License (see VARARGIN)

% Choose default command line output for MARIE_License
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MARIE_License wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = MARIE_License_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


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


% --- Executes during object creation, after setting all properties.
function GPLLogo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GPLLogo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate GPLLogo
FIG = gcbf;
figure(FIG);

if ispc
    % create the LOGO
    logoim = imread('GPLV3_Logo.png');
    image(logoim);
    axis image;
    set(gca,'XTick', []);
    set(gca,'YTick', []);
    
else
    % create the LOGO
    logoim = imread('GPLV3_Logo.png');
    image(logoim);
    axis image;
    set(gca,'XTick', []);
    set(gca,'YTick', []);

end


% --- Executes on button press in Close_Button.
function Close_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Close_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcbf);
close(handles.figure1);
