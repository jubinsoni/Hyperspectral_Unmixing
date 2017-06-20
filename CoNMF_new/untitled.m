function varargout = untitled(varargin)
%UNTITLED M-file for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to untitled_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      UNTITLED('CALLBACK') and UNTITLED('CALLBACK',hObject,...) call the
%      local function named CALLBACK in UNTITLED.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 07-Jun-2017 19:08:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
clc
% Choose default command line output for untitled
handles.output = hObject;
% axes(handles.axes1)
% AA=imread('IIIT.JPG');
% imshow(AA);
% axes(handles.axes2)
% AB=imread('lsi.JPG');
% imshow(AB);
axes(handles.axes4)
AC=imread('empty.png');
imshow(AC);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
movegui('northwest')


% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in uipanel4.
function uipanel4_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel4 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fn1,pn]=uigetfile('*.hdr','select hdr file')
% complete=strcat(pn,fn)
% set(handles.edit1,'string',complete);


i=1;
pn1='';
for i=1:1:size(fn1,2)
    if(fn1(i) == '.')
        break;
    end
     pn1=strcat(pn1,fn1(1,i));
end
setappdata(0,'pn_1',pn1);
setappdata(0,'fn_1',fn1)

% --- Executes during object creation, after setting all properties.
function radiobutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function radiobutton2_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        fn1=getappdata(0,'fn_1')
        pn1=getappdata(0,'pn_1')
        info = envihdrread(fn1);
        D=envidataread(pn1,info);

  % estimated number of  endmembers
  
r800=((D(:,:,44)+D(:,:,45))/2);
r445=((D(:,:,9)+D(:,:,10))/2);
r680=D(:,:,33);
r531=((D(:,:,18)+D(:,:,19))/2);
r570=D(:,:,22);
r415=D(:,:,7);
r435=D(:,:,9);
r700=D(:,:,35);
r695=((D(:,:,34)+D(:,:,35))/2);
r740=(D(:,:,39));
r720=(D(:,:,37));
r750=(D(:,:,40));
r760=(D(:,:,41));
r520=D(:,:,17);
r550=D(:,:,20);
r705=((D(:,:,35)+D(:,:,36))/2);
r717=((D(:,:,36)+D(:,:,37))/2);
r491=((D(:,:,14)+D(:,:,15))/2);
r924=(D(:,:,78));
r703=D(:,:,35);
r423=D(:,:,8);
r693=D(:,:,34);
r1770=D(:,:,162);

rad_on=get(handles.uipanel4,'SelectedObject');
if(rad_on == handles.radiobutton20)
    set(handles.salida,'String','1')
    x1=r800-r445;
    y1=r800-r680;
    var1=x1./y1;
    imagesc(var1)
    set(handles.salida,'String','1')
    dlmwrite('MyFile.txt',var1,'delimiter','\t','precision',3)
elseif(rad_on == handles.radiobutton14)
    x2=r531-r570;
    y2=r531+r570;
    var2=x2./y2;
    imagesc(var2)
    dlmwrite('MyFile.txt',var2,'delimiter','\t','precision',3)
    set(handles.salida,'String','2')
elseif(rad_on == handles.radiobutton15)
    x3=r415-r435;
    y3=r415+r435;
    var3=x3./y3;
    imagesc(var3)
    dlmwrite('MyFile.txt',var3,'delimiter','\t','precision',3)
    set(handles.salida,'String','3')
elseif(rad_on == handles.radiobutton16)
    men=menu('choose a suitable method','Carter index(CI)','Gibelson and merzylah index (GMI)','Vogelman index (VOG)');
    if(men == 1)
        x41=r760;
        y41=r695;
        var41=x41./y41;
        imagesc(var41)
        dlmwrite('MyFile.txt',var41,'delimiter','\t','precision',3)
        set(handles.salida,'String','41')
    elseif(men == 2)
        x42=r750;
        y42=r700;
        var42=x42./y42;
        imagesc(var42)
        dlmwrite('MyFile.txt',var42,'delimiter','\t','precision',3)
        set(handles.salida,'String','42')
    else
        x43=r740;
        y43=r720;
        var43=x43./y43;
        imagesc(var43)
        dlmwrite('MyFile.txt',var43,'delimiter','\t','precision',3)
        set(handles.salida,'String','43')
    end
elseif(rad_on == handles.radiobutton17)
        x5=((r800./r520)-(r800./r550));
        y5=1;
        var5=x5./y5;
        imagesc(var5)
        dlmwrite('MyFile.txt',var5,'delimiter','\t','precision',3)
        set(handles.salida,'String','5')
else
    men=menu('Choose a method','NI_Tian','NI_Wang','NI_Ferwerda')
    if(men == 1)
        x61=r705;
        y61=r717+r491;
        var61=x61./y61;
        imagesc(var61)
        dlmwrite('MyFile.txt',var61,'delimiter','\t','precision',3)
        set(handles.salida,'String','61')
    elseif(men == 2)
        x62=r924-r703+2*(r423);
        y62=r924+r703-2*(r423);
        var62=x62./y62;
        imagesc(var62)
        dlmwrite('MyFile.txt',var62,'delimiter','\t','precision',3)
        set(handles.salida,'String','62')
    else
        x63=r693-r1770;
        y63=r693+r1770;
        var63=x63./y63;
        imagesc(var63)
        dlmwrite('MyFile.txt',var63,'delimiter','\t','precision',3)
        set(handles.salida,'String','63')
    end
end

% --- Executes on button press in radiobutton13.
function radiobutton13_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton13


% --- Executes on button press in radiobutton14.
function radiobutton14_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton14


% --- Executes on button press in radiobutton18.
function radiobutton18_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton18
