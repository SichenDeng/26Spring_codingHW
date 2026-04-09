function varargout = HeatExchangerGUI(varargin)
%HEATEXCHANGERGUI M-file for HeatExchangerGUI.fig
%      HEATEXCHANGERGUI, by itself, creates a new HEATEXCHANGERGUI or raises the existing
%      singleton*.
%
%      H = HEATEXCHANGERGUI returns the handle to a new HEATEXCHANGERGUI or the handle to
%      the existing singleton*.
%
%      HEATEXCHANGERGUI('Property','Value',...) creates a new HEATEXCHANGERGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to HeatExchangerGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      HEATEXCHANGERGUI('CALLBACK') and HEATEXCHANGERGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in HEATEXCHANGERGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HeatExchangerGUI

% Last Modified by GUIDE v2.5 03-Jul-2015 01:42:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HeatExchangerGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @HeatExchangerGUI_OutputFcn, ...
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


% --- Executes just before HeatExchangerGUI is made visible.
function HeatExchangerGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for HeatExchangerGUI
handles.output = hObject;

%initialize image
axes(handles.axes2)
imshow('Parallel.png')

%set variable names with subscripts
[handles.set1, handles.pos1] = javacomponent('javax.swing.JLabel');
set(handles.set1, 'text', '<HTML>m<SUB>h  </SUB>     (<SUP>kg</SUP>/<SUB>s</SUB>)</HTML>')
set(handles.pos1, 'position', [35,500,60,25])

[handles.dot1, handles.posdot1] = javacomponent('javax.swing.JLabel');
set(handles.dot1, 'text', '<HTML>.</HTML>');
set(handles.posdot1, 'position', [35,510,60,25])

[handles.set2, handles.pos2] = javacomponent('javax.swing.JLabel');
set(handles.set2, 'text', '<HTML>m<SUB>c  </SUB>     (<SUP>kg</SUP>/<SUB>s</SUB>)</HTML>')
set(handles.pos2, 'position', [35,464,60,25])

[handles.dot2, handles.posdot2] = javacomponent('javax.swing.JLabel');
set(handles.dot2, 'text', '<HTML>.</HTML>');

[handles.set3, handles.pos3] = javacomponent('javax.swing.JLabel');
set(handles.set3, 'text', '<HTML>c<SUB>p,h  </SUB>    ( <SUP>J</SUP>/<SUB>kg&#8226K</SUB>)</HTML>')
set(handles.pos3, 'position', [35,427,80,25])

[handles.set4, handles.pos4] = javacomponent('javax.swing.JLabel');
set(handles.set4, 'text', '<HTML>c<SUB>p,c  </SUB>    ( <SUP>J</SUP>/<SUB>kg&#8226K</SUB>)</HTML>')
set(handles.pos4, 'position', [35,392,80,25])

[handles.set5, handles.pos5] = javacomponent('javax.swing.JLabel');
set(handles.set5, 'text', '<HTML>T<SUB>h,i  </SUB>      (&#176C)</HTML>')
set(handles.pos5, 'position', [35,355,60,25])

[handles.set6, handles.pos6] = javacomponent('javax.swing.JLabel');
set(handles.set6, 'text', '<HTML>T<SUB>c,i  </SUB>      (&#176C)</HTML>')
set(handles.pos6, 'position', [35,319,60,25])

[handles.set7, handles.pos7] = javacomponent('javax.swing.JLabel');
set(handles.set7, 'text', '<HTML>U   (<SUP>W</SUP>/<SUB>m<SUP>2</SUP>&#8226K</SUB>)</HTML>')
set(handles.pos7, 'position', [35,283,80,25])

[handles.set8, handles.pos8] = javacomponent('javax.swing.JLabel');
set(handles.set8, 'text', '<HTML>L   (m)</HTML>')
set(handles.pos8, 'position', [35,246,50,25])

[handles.set9, handles.pos9] = javacomponent('javax.swing.JLabel');
set(handles.set9, 'text', '<HTML>D   (cm)</HTML>')
set(handles.pos9, 'position', [35,210,50,25])

[handles.set10, handles.pos10] = javacomponent('javax.swing.JLabel');
set(handles.set10, 'text', '<HTML>T<SUB>h,o  </SUB> = </HTML>')
set(handles.pos10, 'position', [560,513,40,30])

[handles.set11, handles.pos11] = javacomponent('javax.swing.JLabel');
set(handles.set11, 'text', '<HTML>T<SUB>c,o  </SUB> = </HTML>')
set(handles.pos11, 'position', [560,478,40,30])

[handles.set12, handles.pos12] = javacomponent('javax.swing.JLabel');
set(handles.set12, 'text', '<HTML>&#949 = </HTML>')
set(handles.pos12, 'position', [560,460,40,30])

figure1_SizeChangedFcn(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HeatExchangerGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = HeatExchangerGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function plot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate plot


% --- Executes during object deletion, before destroying properties.
function plot_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in runbutton.
function runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

mhdot = str2double(get(handles.edit1,'String'));
mcdot = str2double(get(handles.edit2,'String'));
ch = str2double(get(handles.edit3,'String'));
cc = str2double(get(handles.edit4,'String'));
Thi = str2double(get(handles.edit5,'String'));
Tci = str2double(get(handles.edit6,'String'));
U = str2double(get(handles.edit7,'String'));
L = str2double(get(handles.edit8,'String'));
D = str2double(get(handles.edit9,'String'))/100;

if (Tci > Thi)
    errordlg('The hot inlet temperature must be greater than the cold inlet temperature','Error');
    return
end
 
axes(handles.plot)
if get(handles.parallelbutton, 'Value')  
    [Tho, Tco, e] = ParallelFlowFunction(mhdot,mcdot,ch,cc,Thi,Tci,U,L,D);
    
    set(handles.Th, 'String', sprintf('%.2f %cC',Tho, char(176)))
    set(handles.Tc, 'String', sprintf('%.2f %cC',Tco, char(176)))
    set(handles.epsilon, 'String', sprintf('%.2f', e))
else
    [Tho, Tco, e] = CounterFlowFunction(mhdot,mcdot,ch,cc,Thi,Tci,U,L,D);
    
    set(handles.Th, 'String', sprintf('%.2f %cC',Tho, char(176)))
    set(handles.Tc, 'String', sprintf('%.2f %cC',Tco, char(176)))
    set(handles.epsilon, 'String', sprintf('%.2f', e))
end

guidata(hObject,handles)

% --- Executes on button press in parallelbutton.
function parallelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to parallelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parallelbutton
axes(handles.axes2)
imshow('Parallel.png')

guidata(hObject,handles)

% --- Executes on button press in counterbutton.
function counterbutton_Callback(hObject, eventdata, handles)
% hObject    handle to counterbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of counterbutton
axes(handles.axes2)
imshow('Counter.png')

guidata(hObject,handles)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get new figure size
if (handles.figure1.Position(4) < 600)
    handles.figure1.Position(4) = 600;
end
if (handles.figure1.Position(3) < 740)
    handles.figure1.Position(3) = 740;
end

figheight = handles.figure1.Position(4);
figwidth = handles.figure1.Position(3);

%set new positions
%panel
panelbottom = (.15*figheight);
panelleft = .05*figwidth;
panelheight = .8*figheight;
panelwidth = .25*figwidth;
handles.uipanel1.Position = [panelleft panelbottom panelwidth panelheight];

%   edit text fields
editleft = .55*panelwidth;
editwidth = .4*panelwidth;
editheight = .08*panelheight;
editfontsize = .4*editheight;

edit1bottom = .892*panelheight;
handles.edit1.Position = [editleft edit1bottom editwidth editheight];
handles.edit1.FontSize = editfontsize;

edit2bottom = .784*panelheight;
handles.edit2.Position = [editleft edit2bottom editwidth editheight];
handles.edit2.FontSize = editfontsize;

edit3bottom = .676*panelheight;
handles.edit3.Position = [editleft edit3bottom editwidth editheight];
handles.edit3.FontSize = editfontsize;

edit4bottom = .568*panelheight;
handles.edit4.Position = [editleft edit4bottom editwidth editheight];
handles.edit4.FontSize = editfontsize;

edit5bottom = .460*panelheight;
handles.edit5.Position = [editleft edit5bottom editwidth editheight];
handles.edit5.FontSize = editfontsize;

edit6bottom = .352*panelheight;
handles.edit6.Position = [editleft edit6bottom editwidth editheight];
handles.edit6.FontSize = editfontsize;

edit7bottom = .244*panelheight;
handles.edit7.Position = [editleft edit7bottom editwidth editheight];
handles.edit7.FontSize = editfontsize;

edit8bottom = .136*panelheight;
handles.edit8.Position = [editleft edit8bottom editwidth editheight];
handles.edit8.FontSize = editfontsize;

edit9bottom = .028*panelheight;
handles.edit9.Position = [editleft edit9bottom editwidth editheight];
handles.edit9.FontSize = editfontsize;

%   variable names
posleft = panelleft + .05*panelwidth;
poswidth = .4*panelwidth;
posheight = .08*panelheight;

pos1bottom = .892*panelheight + panelbottom;
handles.pos1.Position = [posleft pos1bottom poswidth posheight];

dot1bottom = pos1bottom + .6*posheight;
dot1width = 10;
dot1height = 10;
dot1left = posleft + .45*dot1width;
handles.posdot1.Position = [dot1left dot1bottom dot1width dot1height];

pos2bottom = .784*panelheight + panelbottom;
handles.pos2.Position = [posleft pos2bottom poswidth posheight];

dot2bottom = pos2bottom + .6*posheight;
dot2width = 10;
dot2height = 10;
dot2left = posleft + .45*dot1width;
handles.posdot2.Position = [dot2left dot2bottom dot2width dot2height];

pos3bottom = .676*panelheight + panelbottom;
handles.pos3.Position = [posleft pos3bottom poswidth posheight];

pos4bottom = .568*panelheight + panelbottom;
handles.pos4.Position = [posleft pos4bottom poswidth posheight];

pos5bottom = .460*panelheight + panelbottom;
handles.pos5.Position = [posleft pos5bottom poswidth posheight];

pos6bottom = .352*panelheight + panelbottom;
handles.pos6.Position = [posleft pos6bottom poswidth posheight];

pos7bottom = .244*panelheight + panelbottom;
handles.pos7.Position = [posleft pos7bottom poswidth posheight];

pos8bottom = .136*panelheight + panelbottom;
handles.pos8.Position = [posleft pos8bottom poswidth posheight];

pos9bottom = .028*panelheight + panelbottom;
handles.pos9.Position = [posleft pos9bottom poswidth posheight];

%button
buttonwidth = .5*panelwidth;
buttonheight = panelbottom/2;
buttonleft = panelleft + panelwidth/2 - buttonwidth/2;
buttonbottom = panelbottom/4;
handles.runbutton.Position = [buttonleft buttonbottom buttonwidth buttonheight];

%axes
plotbottom = buttonbottom + buttonheight;
plotleft = panelleft + panelwidth + .08*figwidth;
plotwidth = .57*figwidth;
plotheight = .6*figheight;
handles.plot.Position = [plotleft plotbottom plotwidth plotheight];

%buttongroup
groupwidth = (.8/3)*plotwidth;
groupbottom = .8*figheight;
groupleft = plotleft;
groupheight = .15*figheight;
handles.uibuttongroup1.Position = [groupleft groupbottom groupwidth groupheight];

parallelwidth = groupwidth;
parallelleft = 0;
parallelbottom = .55*groupheight;
parallelheight = .35*groupheight;
handles.parallelbutton.Position = [parallelleft parallelbottom parallelwidth parallelheight];
handles.parallelbutton.FontSize = .4*parallelheight;

counterwidth = groupwidth;
counterleft = 0;
counterbottom = .1*groupheight;
counterheight = .35*groupheight;
handles.counterbutton.Position = [counterleft counterbottom counterwidth counterheight];
handles.counterbutton.FontSize = .4*counterheight;

%picture
picbottom = groupbottom;
picleft = groupleft + groupwidth + .1*plotwidth;
picwidth = (.8/3)*plotwidth;
picheight = .15*figheight;
handles.axes2.Position = [picleft picbottom picwidth picheight];

%temperature out
pos10bottom = .92*figheight;
pos10left = picleft + picwidth + .1*plotwidth;
pos10width = .1*plotwidth;
pos10height = .03*figheight;
handles.pos10.Position = [pos10left pos10bottom pos10width pos10height];

pos11height = .03*figheight;
pos11bottom = .86*figheight;
pos11left = pos10left;
pos11width = pos10width;
handles.pos11.Position = [pos11left pos11bottom pos11width pos11height];

thleft = pos10left + pos10width;
thbottom = pos10bottom;
thwidth = plotleft + plotwidth - pos10left - pos10width;
thheight = pos10height;
handles.Th.Position = [thleft thbottom thwidth thheight];

tcleft = thleft;
tcbottom = pos11bottom;
tcwidth = thwidth;
tcheight = thheight;
handles.Tc.Position = [tcleft tcbottom tcwidth tcheight];

%effectiveness
pos12left = pos10left;
pos12bottom = .8*figheight;
pos12width = pos10width;
pos12height = pos10height;
handles.pos12.Position = [pos12left pos12bottom pos12width pos12height];

epsilonwidth = thwidth;
epsilonleft = thleft;
epsilonbottom = pos12bottom;
epsilonheight = thheight;
handles.epsilon.Position = [epsilonleft epsilonbottom epsilonwidth epsilonheight];
