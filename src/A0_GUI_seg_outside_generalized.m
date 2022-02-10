%guide A0_GUI_seg_outside_generalized.fig; use to open GUI figure to edit
function varargout = A0_GUI_seg_outside_generalized(varargin)
% A0_GUI_SEG_OUTSIDE_GENERALIZED MATLAB code for A0_GUI_seg_outside_generalized.fig
%      A0_GUI_SEG_OUTSIDE_GENERALIZED, by itself, creates a new A0_GUI_SEG_OUTSIDE_GENERALIZED or raises the existing
%      singleton*.
%
%      H = A0_GUI_SEG_OUTSIDE_GENERALIZED returns the handle to a new A0_GUI_SEG_OUTSIDE_GENERALIZED or the handle to
%      the existing singleton*.
%
%      A0_GUI_SEG_OUTSIDE_GENERALIZED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in A0_GUI_SEG_OUTSIDE_GENERALIZED.M with the given input arguments.
%
%      A0_GUI_SEG_OUTSIDE_GENERALIZED('Property','Value',...) creates a new A0_GUI_SEG_OUTSIDE_GENERALIZED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before A0_GUI_seg_outside_generalized_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.     All inputs are passed to A0_GUI_seg_outside_generalized_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help A0_GUI_seg_outside_generalized

% Last Modified by GUIDE v2.5 03-Feb-2022 13:51:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @A0_GUI_seg_outside_generalized_OpeningFcn, ...
                   'gui_OutputFcn',  @A0_GUI_seg_outside_generalized_OutputFcn, ...
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


% --- Executes just before A0_GUI_seg_outside_generalized is made visible.
function A0_GUI_seg_outside_generalized_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to A0_GUI_seg_outside_generalized (see VARARGIN)
% Choose default command line output for A0_GUI_seg_outside_generalized
handles.output = hObject;
%Below are defaults for every variable
handles.thA = [1,1,200]; % Thresholds for RNA spot intensities
handles.ths  = [50 50 21 50 101 50 50 50 50]; % Thresholds for RNA spot intensities % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS defines which fluorecnet dyes have been imaged 
handles.Ych = [0 0 0 0 0 0 0 0 0] ; %first is the total number of channels, second is the marker slice, and third is the boundary slice 
handles.img_stacks = [0,0];
handles.Yth = 0;
handles.outfile_prefix_RNA = '';
handles.ManTh = 0;
handles.Ywin = 0;
handles.Yim = 0;
handles.file_type = 'Multiple';
handles.segment_mode = 'last5';
handles.Seg_thres = 0;
handles.fish_dir = '';
handles.outfile_prefix_seg = ''; %Where the siles will be saved to
handles.max_int_thres = 0;
handles.max_int_spot_th = 0;
handles.min_nucleus_size = 0;
handles.max_nucleus_size = 0;
handles.min_cell_size = 0;
handles.max_cell_size = 0;
handles.yeast_seg = 0;
handles.images1 = {};
handles.im_prefixes = {};
handles.input_name = '';
handles.node_num = 1;
handles.Diffstack = 0;
handles.DAPI_slice = 1;
handles.TRANS_slice = 1;
handles.DAPI_min= 0;
handles.DAPI_max = 1;
handles.TRANS_min = 0;
handles.TRANS_max = 1;
handles.Xlim_fig = [0 2048];
handles.Ylim_fig = [0 2048];
handles.nuc_border = 0;
handles.cell_border = 0;
handles.thicken_num = 0;
handles.show_nuc_bord = 1;
handles.show_cell_bord = 1;
handles.DAPI_prefix = '';
handles.TRANS_prefix = '';
handles.images_DAPI = {};
handles.images_TRANS = {};
handles.Diffstack = 0;
handles.Im_name = 'Sample Image';
set(handles.Min_Nuc_edit,'Tooltip',...
strcat('Type the number of pixels in what you think will be ', char(10),...
'the largest marker size (ex: largest nucleus size in', char(10), ...
'DAPI). Try to be exact or a little larger', char(10),... 
'You can use the section to the right to load and', char(10), char(10), ...
'measure, or you can use another program such as ImageJ.'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes A0_GUI_seg_outside_generalized wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = A0_GUI_seg_outside_generalized_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Im_Dir_Push.
function Im_Dir_Push_Callback(hObject, eventdata, handles)
% hObject    handle to Im_Dir_Push (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.fish_dir= uigetdir();
Ywin = handles.Ywin;
if Ywin
    handles.fish_dir = [handles.fish_dir '\'];
else
    handles.fish_dir = [handles.fish_dir '/'];
end
fish_dir = handles.fish_dir;
set(handles.Im_Dir_edit,'String',fish_dir);
% if handles.Diffstack
%  [images1,folders1] = Get_Dir_All(handles.fish_dir,handles.Ywin);
%  clear images2   
%  images2 = {};
%  for i = 1:size(images1,2)
%      if Ywin
%      beginning1 = strfind(images1{i},'\');
%      else
%          beginning1 = strfind(images1{i},'/');
%      end
%      beginning1 = beginning1(end)+1;    %Find the part that is not directory
%      if strfind(images1{i},'_MMStack')
%          end1 = strfind(images1{i},'_MMStack');
%      else
%          end1 = strfind(images1{i},'.');
%      end
%      end1 = end1(1)-1;
%      images2{i} = images1{i}(beginning1:end1);
%  end 
%   handles.images1 = images1;
%   handles.images2 = images2;
%   counter_DAPI = 1;
%   counter_TRANS = 1; 
%   handles.DAPI_prefix
%   for i = 1:size(images2,2)
%       if strfind(images2{i},handles.DAPI_prefix)
%           handles.images_DAPI{counter_DAPI} = images1{i};
%           counter_DAPI = counter_DAPI+1;
%       elseif strfind(images2{i},handles.TRANS_prefix)
%           handles.images_TRANS{counter_TRANS} = images1{i};
%           counter_TRANS = counter_TRANS+1;
%       end
%   end 
% else  
% [images1,im_prefixes] = Get_Dir_ImSt(handles.fish_dir,handles.Ywin);                                           %Finds the file names within the fish directory
%     handles.images1 = images1;
%     handles.im_prefixes = im_prefixes;
% end
[images1,folders1] = Get_Dir_All(handles.fish_dir,handles.Ywin);
 clear images2 im_prefixes  
 images2 = {};          %This will have the file names (not the entire directory and not the type of file)
 for i = 1:size(images1,2)
     if Ywin                     %The type of slash depends on whether it is windows
     beginning1 = strfind(images1{i},'\');      
     else
         beginning1 = strfind(images1{i},'/');
     end
     beginning1 = beginning1(end)+1;    %Find the part that is not directory
     if strfind(images1{i},'_MMStack')
         end1 = strfind(images1{i},'_MMStack'); %look for MMStack since it is often at the end of the file names
     elseif strfind(images1{i},'_MMImages')
         end1 = strfind(images1{i},'_MMImages'); %look for MMImages since it is often at the end of the file names
     elseif strfind(images1{i},'.ome')
         end1 = strfind(images1{i},'.ome');        %Look for a period otherwise as the end of the file name
     else
         end1 = strfind(images1{i},'.');        %Look for a period otherwise as the end of the file name
     end
     end1 = end1(1)-1;
     images2{i} = images1{i}(beginning1:end1);  %Populate the cell array with the file names
 end 
 images2';  %For visualizing all files extracted (remove semicolon)
  handles.images1 = images1;
  handles.images2 = images2;
if handles.Diffstack
    clear handles.images_DAPI handles.images_TRANS
    handles.images_DAPI = {};
    handles.images_DAPI = {};
  counter_DAPI = 1;             %counter for marker images found +1 (used to populate both DAPI and TRANS cell arrays)
  counter_TRANS = 0;            %counter for TRANS images found (used to determine if TRANS prefix was correct)
  handles.DAPI_prefix
  for i = 1:size(images2,2) 
      if strfind(images2{i},handles.DAPI_prefix) == 1
          handles.images_DAPI{counter_DAPI} = images1{i};
          im_prefixes{counter_DAPI} =  images2{i}(size(handles.DAPI_prefix,2)+1:end);   %populate imprefixes (used for later segmentation file names)
          for j = 1:size(images2,2)
              if strfind(images2{j},handles.TRANS_prefix) == 1 &...         %Look for TRANS prefix at beginning
                      strcmp(images2{j}(size(handles.TRANS_prefix,2)+1:end),... %Look for the rest of image name to match between TRANS
                      images2{i}(size(handles.DAPI_prefix,2)+1:end))        %and DAPI
                      handles.images_TRANS{counter_DAPI} = images1{j}; 
              end
          end
          counter_TRANS = counter_TRANS+1;                                  %Add to counter (even if match not found)
          counter_DAPI = counter_DAPI+1;                                    %Add to counter (even if match not found)
       if strfind(images2{i},handles.TRANS_prefix) %searching for TRANS prefix
%            handles.images_TRANS{counter_TRANS} = images1{i};
           counter_TRANS = counter_TRANS+1;
       end
      end
  end
  if counter_TRANS == 0
      msgbox('There were no images found with the boundary prefix')
  end
  if isempty(handles.images_DAPI)
      msgbox('There were no images found with the marker prefix')     
  elseif isempty(handles.images_TRANS)
      msgbox('There were no images paired with the prefixes specified')
  end
else
%[images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                %Finds the file names within the fish directory (old way)
im_prefixes = images2;
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
end
images1
images2
folders1
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Yim_checkbox.
function Yim_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Yim_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles = guidata(hObject);
handles.Yim = get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of Yim_checkbox
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Ywin_checkbox.
function Ywin_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Ywin_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ywin = get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of Ywin_checkbox
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in File_Type_popupmenu.
% function File_Type_popupmenu_Callback(hObject, eventdata, handles);
% % hObject    handle to File_Type_popupmenu (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: contents = cellstr(get(hObject,'String')) returns File_Type_popupmenu contents as cell array
% %        contents{get(hObject,'Value')} returns selected item from File_Type_popupmenu
% handles.contents = cellstr(get(hObject,'String')); %returns File_Type_popupmenu contents as cell array
% handles.file_type = handles.contents{get(hObject,'Value')};
% % Update handles structure
% guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
% function File_Type_popupmenu_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to File_Type_popupmenu (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% % Hint: popupmenu controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes on button press in Max_Int_Spot_Thres_checkbox.
function Max_Int_Spot_Thres_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Max_Int_Spot_Thres_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.max_int_thres = get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of Max_Int_Spot_Thres_checkbox
% Update handles structure
guidata(hObject, handles);



function Min_Nuc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Min_Nuc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.min_nucleus_size = str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Min_Nuc_edit as text
%        str2double(get(hObject,'String')) returns contents of Min_Nuc_edit as a double


% --- Executes during object creation, after setting all properties.
function Min_Nuc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_Nuc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Max_Nuc_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Max_Nuc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.max_nucleus_size = str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of Max_Nuc_edit as text
%        str2double(get(hObject,'String')) returns contents of Max_Nuc_edit as a double


% --- Executes during object creation, after setting all properties.
function Max_Nuc_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_Nuc_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Min_Cell_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Min_Cell_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.min_cell_size = str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of Min_Cell_edit as text
%        str2double(get(hObject,'String')) returns contents of Min_Cell_edit as a double


% --- Executes during object creation, after setting all properties.
function Min_Cell_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_Cell_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Max_Cell_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Max_Cell_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.max_cell_size = str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of Max_Cell_edit as text
%        str2double(get(hObject,'String')) returns contents of Max_Cell_edit as a double


% --- Executes during object creation, after setting all properties.
function Max_Cell_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_Cell_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Seg_Dir_pushbutton.
function Seg_Dir_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Seg_Dir_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.outfile_prefix_seg= uigetdir();
if handles.Ywin
handles.outfile_prefix_seg = [handles.outfile_prefix_seg '\'];
else
    handles.outfile_prefix_seg = [handles.outfile_prefix_seg '/'];
end
set(handles.Seg_Dir_edit,'String',handles.outfile_prefix_seg);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in Save_In_pushbutton.
function Save_In_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Save_In_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;
ths = handles.ths;
outfile_prefix_RNA = handles.outfile_prefix_RNA;
Yth = handles.Yth;
Ych = handles.Ych;
ManTh = handles.ManTh;
Ywin = handles.Ywin;
Yim = handles.Yim;
file_type = handles.file_type;
segment_mode = handles.segment_mode;
Seg_thres = handles.Seg_thres;
fish_dir = handles.fish_dir;
outfile_prefix_seg = handles.outfile_prefix_seg;
max_int_thres = handles.max_int_thres;
max_int_spot_th = handles.max_int_spot_th;
min_nucleus_size = handles.min_nucleus_size;
max_nucleus_size = handles.max_nucleus_size;
min_cell_size = handles.min_cell_size;
max_cell_size = handles.max_cell_size;
img_stacks = handles.img_stacks;
yeast_seg = handles.yeast_seg;
input_name = handles.input_name;
Diffstack = handles.Diffstack;
DAPI_prefix = handles.DAPI_prefix;
TRANS_prefix = handles.TRANS_prefix; 
save(input_name,'outfile_prefix_RNA','thA','ths','Yth','Ych','ManTh','Ywin','Yim',...
'file_type','segment_mode','Seg_thres','fish_dir','outfile_prefix_seg','min_nucleus_size',...
'max_nucleus_size','min_cell_size','max_cell_size','img_stacks','yeast_seg','Diffstack','DAPI_prefix','TRANS_prefix');

% --- Executes on button press in Load_In_pushbutton.
function Load_In_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Load_In_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Loads variables from file and gives them to handles
load(get(handles.Load_In_edit,'String'));
try handles.thA = thA; end % Thresholds for RNA spot intensities 
try handles.ths  = ths; end % Thresholds for RNA spot intensities % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS defines which fluorecnet dyes have been imaged 
try handles.Ych = Ych ; end
try handles.img_stacks = img_stacks; end
try handles.Yth = Yth; end
% try handles.outfile_prefix_RNA = outfile_prefix_RNA; end
try handles.ManTh = ManTh; end
try handles.Ywin = Ywin; end
try handles.Yim = Yim; end
try handles.file_type = file_type; end
try handles.segment_mode = segment_mode; end
try handles.Seg_thres = Seg_thres; end
try handles.fish_dir = fish_dir; end
try handles.outfile_prefix_seg = outfile_prefix_seg; end
try handles.max_int_thres = max_int_thres; end
try handles.max_int_spot_th = max_int_spot_th; end
try handles.min_nucleus_size = min_nucleus_size; end
try handles.max_nucleus_size = max_nucleus_size; end
try handles.min_cell_size = min_cell_size; end
try handles.max_cell_size = max_cell_size; end
try handles.yeast_seg = yeast_seg; end
try handles.input_name = get(handles.Load_In_edit,'String'); end
try handles.Diffstack = Diffstack; catch Diffstack = 0; handles.Diffstack = Diffstack; end
try handles.TRANS_prefix = TRANS_prefix; catch TRANS_prefix = ''; handles.TRANS_prefix = TRANS_prefix; end
try handles.DAPI_prefix = DAPI_prefix; catch DAPI_prefix = ''; handles.DAPI_prefix = TRANS_prefix; end
%% Populate everything on the form with updated variables
% set(handles.Beg_thres_edit,'Value',handles.thA(1)); % Thresholds for RNA spot intensities
% set(handles.Int_edit,'Value',handles.thA(2));
% set(handles.End_Thres_edit,'Value',handles.thA(3));
%set(handles.handles.ths  % Thresholds for RNA spot intensities % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS defines which fluorecnet dyes have been imaged 
% set(handles.DAPI_checkbox,'Value',handles.Ych(2)) ;
% set(handles.CY5_checkbox,'Value',handles.Ych(3));
% set(handles.AF594_checkbox,'Value',handles.Ych(4)) ;
% set(handles.TMR_checkbox,'Value',handles.Ych(5)) ;
% set(handles.YFP_checkbox,'Value',handles.Ych(6)) ;
% set(handles.GFP_checkbox,'Value',handles.Ych(7)) ;
% set(handles.AF700_checkbox,'Value',handles.Ych(8)); 
% set(handles.TRANS_checkbox,'Value',handles.Ych(9)) ;
set(handles.Diff_Stack_checkbox,'Value',Diffstack)
%%%Make changes depending on Diffstack (whether DAPI and TRANS images are
%%%in separate stacks)
if handles.Diffstack == 0
    set(handles.Mark_im_prefix_edit, 'enable', 'off')
    set(handles.Bound_im_prefix_edit, 'enable', 'off')
    set(handles.text18, 'enable', 'off')
    set(handles.text17, 'enable', 'off')  
    set(handles.Load_Sample_pushbutton, 'enable', 'on')
    set(handles.Boundary_slice_edit, 'enable', 'on') 
    set(handles.text19, 'enable', 'on') 
    set(handles.Marker_slice_edit, 'enable', 'on') 
    set(handles.text20, 'enable', 'on') 
    set(handles.Total_channel_edit, 'enable', 'on') 
    set(handles.text21, 'enable', 'on')       
elseif handles.Diffstack == 1
    set(handles.Mark_im_prefix_edit, 'enable', 'on')
    set(handles.Bound_im_prefix_edit, 'enable', 'on')
    set(handles.text18, 'enable', 'on')
    set(handles.text17, 'enable', 'on')
    set(handles.Load_Sample_pushbutton, 'enable', 'off')
    set(handles.Boundary_slice_edit, 'enable', 'off') 
    set(handles.text19, 'enable', 'off') 
    set(handles.Marker_slice_edit, 'enable', 'off') 
    set(handles.text20, 'enable', 'off') 
    set(handles.Total_channel_edit, 'enable', 'off') 
    set(handles.text21, 'enable', 'off')      
end
set(handles.Mark_im_prefix_edit,'String',DAPI_prefix)
set(handles.Bound_im_prefix_edit,'String',TRANS_prefix)
set(handles.Low_Slice_edit,'String',num2str(handles.img_stacks(1)));
set(handles.High_Slice_edit,'String',num2str(handles.img_stacks(2)));
set(handles.ManTh_checkbox,'Value',handles.ManTh);
%set(handles.handles.Yth
% set(handles.RNA_Dir_edit,'String',handles.outfile_prefix_RNA);
%set(handles.handles.ManTh
set(handles.Ywin_checkbox,'Value',handles.Ywin);
set(handles.Yim_checkbox,'Value',handles.Yim);
% handles.contents = cellstr(get(handles.File_Type_popupmenu,'String')); %returns File_Type_popupmenu contents as cell array
% temp_ind = 1;    %Will have the index of the file type for showing it. but defaults to "File type"
% for i = 1:size(handles.contents)
%     if strcmp(handles.file_type,handles.contents{i})    %compare file type to possible ones listed
%         temp_ind = i;
%     end
% end
% set(handles.File_Type_popupmenu,'Value',temp_ind); %set the file type in the dropdown menu to the correct one  
%set(handles.handles.segment_mode
set(handles.Single_Thres_checkbox,'Value',handles.Seg_thres);
set(handles.Im_Dir_edit,'String',handles.fish_dir);
set(handles.Seg_Dir_edit,'String',handles.outfile_prefix_seg);
% set(handles.Max_Int_Thres_checkbox,'Value',handles.max_int_thres);
% set(handles.Max_Int_Spot_Thres_checkbox,'Value',handles.max_int_spot_th);
set(handles.Min_Nuc_edit,'String',num2str(handles.min_nucleus_size));
set(handles.Max_Nuc_edit,'String',num2str(handles.max_nucleus_size));
set(handles.Min_Cell_edit,'String',num2str(handles.min_cell_size));
set(handles.Max_Cell_edit,'String',num2str(handles.max_cell_size));
set(handles.Total_channel_edit,'String',num2str(handles.Ych(1)));
set(handles.Boundary_slice_edit,'String',num2str(handles.Ych(3)));
set(handles.Marker_slice_edit,'String',num2str(handles.Ych(2)));
%handles.yeast_seg
set(handles.Input_Name_edit,'String',get(handles.Load_In_edit,'String'));
% [images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin) ;   %Gets images if not Diffstack
% handles.images1 = images1;
% handles.im_prefixes = im_prefixes;
% if handles.Diffstack
%  [images1,folders1] = Get_Dir_All(handles.fish_dir,handles.Ywin);
%  clear images2   
%  images2 = {};
%  for i = 1:size(images1,2)
%      if Ywin
%      beginning1 = strfind(images1{i},'\');
%      else
%          beginning1 = strfind(images1{i},'/');
%      end
%      beginning1 = beginning1(end)+1;    %Find the part that is not directory
%      if strfind(images1{i},'_MMStack')
%          end1 = strfind(images1{i},'_MMStack');
%      else
%          end1 = strfind(images1{i},'.');
%      end
%      end1 = end1(1)-1;
%      images2{i} = images1{i}(beginning1:end1);
%  end 
%   handles.images1 = images1;
%   handles.images2 = images2;
%   counter_DAPI = 1;
%   counter_TRANS = 1; 
%   handles.DAPI_prefix
%   for i = 1:size(images2,2)
%       if strfind(images2{i},handles.DAPI_prefix)
%           handles.images_DAPI{counter_DAPI} = images1{i};
%           counter_DAPI = counter_DAPI+1;
%       elseif strfind(images2{i},handles.TRANS_prefix)
%           handles.images_TRANS{counter_TRANS} = images1{i};
%           counter_TRANS = counter_TRANS+1;
%       end
%   end 
% else  
% [images1,im_prefixes] = Get_Dir_ImSt(handles.fish_dir,handles.Ywin);                                           %Finds the file names within the fish directory
%     handles.images1 = images1;
%     handles.im_prefixes = im_prefixes;
% end
[images1,folders1] = Get_Dir_All(handles.fish_dir,handles.Ywin);
 clear images2 im_prefixes  
 images2 = {};          %This will have the file names (not the entire directory and not the type of file)
 for i = 1:size(images1,2)
     if Ywin                     %The type of slash depends on whether it is windows
     beginning1 = strfind(images1{i},'\');      
     else
         beginning1 = strfind(images1{i},'/');
     end
     beginning1 = beginning1(end)+1;    %Find the part that is not directory
     if strfind(images1{i},'_MMStack')
         end1 = strfind(images1{i},'_MMStack'); %look for MMStack since it is often at the end of the file names
     elseif strfind(images1{i},'_MMImages')
         end1 = strfind(images1{i},'_MMImages'); %look for MMImages since it is often at the end of the file names
     elseif strfind(images1{i},'.ome')
         end1 = strfind(images1{i},'.ome');        %Look for a period otherwise as the end of the file name
     else
         end1 = strfind(images1{i},'.');        %Look for a period otherwise as the end of the file name
     end
     end1 = end1(1)-1;
     images2{i} = images1{i}(beginning1:end1);  %Populate the cell array with the file names
 end 
 images2';  %For visualizing all files extracted (remove semicolon)
  handles.images1 = images1;
  handles.images2 = images2;
if handles.Diffstack
    clear handles.images_DAPI handles.images_TRANS
    handles.images_DAPI = {};
    handles.images_DAPI = {};
  counter_DAPI = 1;             %counter for marker images found +1 (used to populate both DAPI and TRANS cell arrays)
  counter_TRANS = 0;            %counter for TRANS images found (used to determine if TRANS prefix was correct)
  handles.DAPI_prefix
  for i = 1:size(images2,2) 
      if strfind(images2{i},handles.DAPI_prefix) == 1
          handles.images_DAPI{counter_DAPI} = images1{i};
          im_prefixes{counter_DAPI} =  images2{i}(size(handles.DAPI_prefix,2)+1:end);   %populate imprefixes (used for later segmentation file names)
          for j = 1:size(images2,2)
              if strfind(images2{j},handles.TRANS_prefix) == 1 &...         %Look for TRANS prefix at beginning
                      strcmp(images2{j}(size(handles.TRANS_prefix,2)+1:end),... %Look for the rest of image name to match between TRANS
                      images2{i}(size(handles.DAPI_prefix,2)+1:end))        %and DAPI
                      handles.images_TRANS{counter_DAPI} = images1{j}; 
              end
          end
          counter_TRANS = counter_TRANS+1;                                  %Add to counter (even if match not found)
          counter_DAPI = counter_DAPI+1;                                    %Add to counter (even if match not found)
       if strfind(images2{i},handles.TRANS_prefix) %searching for TRANS prefix
%            handles.images_TRANS{counter_TRANS} = images1{i};
           counter_TRANS = counter_TRANS+1;
       end
      end
  end
  if counter_TRANS == 0
      msgbox('There were no images found with the boundary prefix')
  end
  if isempty(handles.images_DAPI)
      msgbox('There were no images found with the marker prefix')     
  elseif isempty(handles.images_TRANS)
      msgbox('There were no images paired with the prefixes specified')
  end
else
%[images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                %Finds the file names within the fish directory (old way)
im_prefixes = images2;
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
end
% Update handles structure
guidata(hObject, handles);




% --- Executes on button press in Load_Im_pushbutton.
function Load_Im_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Im_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% [images1,im_prefixes] = Get_Dir_ImSt(handles.fish_dir,handles.Ywin) ;                                          %Finds the file names within the fish directory
% if Ywin == 0
%     handles.images1 = {images1{node_num+1}};
%     handles.im_prefixes = {im_prefixes{node_num+1}};
% end
% Update handles structure
guidata(hObject, handles);


function Im_Dir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Im_Dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Im_Dir_edit as text
%        str2double(get(hObject,'String')) returns contents of Im_Dir_edit as a double
handles.fish_dir = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Im_Dir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Im_Dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Seg_Dir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Seg_Dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.outfile_prefix_seg= uigetdir();
% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Seg_Dir_edit as text
%        str2double(get(hObject,'String')) returns contents of Seg_Dir_edit as a double


% --- Executes during object creation, after setting all properties.
function Seg_Dir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Seg_Dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DAPI_checkbox.
function DAPI_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to DAPI_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ych(2) = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of DAPI_checkbox


% --- Executes on button press in CY5_checkbox.
function CY5_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to CY5_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ych(3) = get(hObject,'Value'); % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of CY5_checkbox


% --- Executes on button press in AF594_checkbox.
function AF594_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to AF594_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ych(4) = get(hObject,'Value'); % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of AF594_checkbox


% --- Executes on button press in TMR_checkbox.
function TMR_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to TMR_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ych(5) = get(hObject,'Value'); % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of TMR_checkbox


% --- Executes on button press in YFP_checkbox.
function YFP_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to YFP_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ych(6) = get(hObject,'Value'); % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of YFP_checkbox


% --- Executes on button press in GFP_checkbox.
function GFP_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to GFP_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ych(7) = get(hObject,'Value'); % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of GFP_checkbox


% --- Executes on button press in AF700_checkbox.
function AF700_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to AF700_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ych(8) = get(hObject,'Value'); % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of AF700_checkbox


% --- Executes on button press in TRANS_checkbox.
function TRANS_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to TRANS_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Ych(9) = get(hObject,'Value'); % YCY7 DAPI CY5 AF594 TMR YFP GFP/AF488 AF700 TRANS; defines which fluorescent dyes have been imaged; 
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of TRANS_checkbox


% --- Executes on button press in Seg_Start_pushbutton.
function Seg_Start_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Seg_Start_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;
ths = handles.ths;
Yth = handles.Yth;
Ych = handles.Ych;
ManTh = handles.ManTh;
% if ManTh
%     Yth = 1;
% end
Ywin = handles.Ywin;
Yim = handles.Yim;
file_type = handles.file_type;
segment_mode = handles.segment_mode;
Seg_thres = handles.Seg_thres;
fish_dir = handles.fish_dir;
outfile_prefix_seg = handles.outfile_prefix_seg;
node_num = handles.node_num;
max_int_thres = handles.max_int_thres;
max_int_spot_th = handles.max_int_spot_th;
min_nucleus_size = handles.min_nucleus_size;
max_nucleus_size = handles.max_nucleus_size;
min_cell_size = handles.min_cell_size;
max_cell_size = handles.max_cell_size;
img_stacks = handles.img_stacks;
yeast_seg = handles.yeast_seg;
images1 = handles.images1;
im_prefixes = handles.im_prefixes;
node_num = handles.node_num;
Diffstack=handles.Diffstack;
fish_dir = handles.fish_dir;
Diffstack
[images1,folders1] = Get_Dir_All(fish_dir,Ywin);
 clear images2 im_prefixes  
 images2 = {};          %This will have the file names (not the entire directory and not the type of file)
 for i = 1:size(images1,2)
     if Ywin                     %The type of slash depends on whether it is windows
     beginning1 = strfind(images1{i},'\');      
     else
         beginning1 = strfind(images1{i},'/');
     end
     beginning1 = beginning1(end)+1;    %Find the part that is not directory
     if strfind(images1{i},'_MMStack')
         end1 = strfind(images1{i},'_MMStack'); %look for MMStack since it is often at the end of the file names
     elseif strfind(images1{i},'_MMImages')
         end1 = strfind(images1{i},'_MMImages'); %look for MMImages since it is often at the end of the file names
     elseif strfind(images1{i},'.ome')
         end1 = strfind(images1{i},'.ome');        %Look for a period otherwise as the end of the file name
     else
         end1 = strfind(images1{i},'.');        %Look for a period otherwise as the end of the file name
     end
     end1 = end1(1)-1;
     images2{i} = images1{i}(beginning1:end1);  %Populate the cell array with the file names
 end 
 images2';  %For visualizing all files extracted (remove semicolon)
  handles.images1 = images1;
  handles.images2 = images2;
if Diffstack
    clear handles.images_DAPI handles.images_TRANS
    handles.images_DAPI = {};
    handles.images_DAPI = {};
  counter_DAPI = 1;             %counter for marker images found +1 (used to populate both DAPI and TRANS cell arrays)
  counter_TRANS = 0;            %counter for TRANS images found (used to determine if TRANS prefix was correct)
  handles.DAPI_prefix
  for i = 1:size(images2,2) 
      if strfind(images2{i},handles.DAPI_prefix) == 1
          handles.images_DAPI{counter_DAPI} = images1{i};
          im_prefixes{counter_DAPI} =  images2{i}(size(handles.DAPI_prefix,2)+1:end);   %populate imprefixes (used for later segmentation file names)
          for j = 1:size(images2,2)
              if strfind(images2{j},handles.TRANS_prefix) == 1 &...         %Look for TRANS prefix at beginning
                      strcmp(images2{j}(size(handles.TRANS_prefix,2)+1:end),... %Look for the rest of image name to match between TRANS
                      images2{i}(size(handles.DAPI_prefix,2)+1:end))        %and DAPI
                      handles.images_TRANS{counter_DAPI} = images1{j}; 
              end
          end
          counter_TRANS = counter_TRANS+1;                                  %Add to counter (even if match not found)
          counter_DAPI = counter_DAPI+1;                                    %Add to counter (even if match not found)
       if strfind(images2{i},handles.TRANS_prefix) %searching for TRANS prefix
%            handles.images_TRANS{counter_TRANS} = images1{i};
           counter_TRANS = counter_TRANS+1;
       end
      end
  end
  if counter_TRANS == 0
      msgbox('There were no images found with the boundary prefix')
  end
  if isempty(handles.images_DAPI)
      msgbox('There were no images found with the marker prefix')     
  elseif isempty(handles.images_TRANS)
      msgbox('There were no images paired with the prefixes specified')
  end
else
%[images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                %Finds the file names within the fish directory (old way)
im_prefixes = images2;
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
end
dxy1 = round(sqrt(max_nucleus_size/pi));  
DAPI_images = handles.images_DAPI
TRANS_images= handles.images_TRANS
A1_segment_predefined_variables_streamlined_generalized(Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_stacks,yeast_seg,images1,im_prefixes,dxy1,Diffstack,DAPI_images,TRANS_images)
msgbox('Segmentation Finished')
% if ManTh
%     Yth = 0;
%     A1_segment_predefined_variables_streamlined_generalized(Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_stacks,yeast_seg,images1,im_prefixes,dxy1,Diffstack,DAPI_images,TRANS_images)
% end
    

% --- Executes on button press in Start_Thresh_pushbutton.
function Start_Thresh_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Start_Thresh_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;
ths = handles.ths;
outfile_prefix_RNA = handles.outfile_prefix_RNA;
Yth = handles.Yth;
Ych = handles.Ych;
ManTh = handles.ManTh;
Ywin = handles.Ywin;
Yim = handles.Yim;
file_type = handles.file_type;
segment_mode = handles.segment_mode;
Seg_thres = handles.Seg_thres;
fish_dir = handles.fish_dir;
outfile_prefix_seg = handles.outfile_prefix_seg;
node_num = handles.node_num;
max_int_thres = handles.max_int_thres;
max_int_spot_th = handles.max_int_spot_th;
Run_RNA_detection_predefined_variables(thA,ths,outfile_prefix_RNA,Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,node_num,max_int_thres,max_int_spot_th)


% --- Executes on button press in Single_Thres_checkbox.
function Single_Thres_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Single_Thres_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Seg_thres = get(hObject,'Value');
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Single_Thres_checkbox



function Input_Name_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Input_Name_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.input_name = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Input_Name_edit as text
%        str2double(get(hObject,'String')) returns contents of Input_Name_edit as a double

% --- Executes during object creation, after setting all properties.
function Input_Name_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Input_Name_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Low_Slice_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Low_Slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
'This ran';
handles.img_stacks(1) =  str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles)
% Hints: get(hObject,'String') returns contents of Low_Slice_edit as text
%        str2double(get(hObject,'String')) returns contents of Low_Slice_edit as a double


% --- Executes during object creation, after setting all properties.
function Low_Slice_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Low_Slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function High_Slice_edit_Callback(hObject, eventdata, handles)
% hObject    handle to High_Slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.img_stacks(2) =  str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of High_Slice_edit as text
%        str2double(get(hObject,'String')) returns contents of High_Slice_edit as a double


% --- Executes during object creation, after setting all properties.
function High_Slice_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to High_Slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Load_In_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Load_In_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Load_In_edit as text
%        str2double(get(hObject,'String')) returns contents of Load_In_edit as a double


% --- Executes during object creation, after setting all properties.
function Load_In_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Load_In_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Load_In_Browse_pushbutton.
function Load_In_Browse_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Load_In_Browse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile;
set(handles.Load_In_edit,'String',[path file]);
load(get(handles.Load_In_edit,'String'));
try handles.thA = thA; end % Thresholds for RNA spot intensities 
try handles.ths  = ths; end % Thresholds for RNA spot intensities % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS defines which fluorecnet dyes have been imaged 
try handles.Ych = Ych ; end
try handles.img_stacks = img_stacks; end
try handles.Yth = Yth; end
% try handles.outfile_prefix_RNA = outfile_prefix_RNA; end
try handles.ManTh = ManTh; end
try handles.Ywin = Ywin; end
try handles.Yim = Yim; end
try handles.file_type = file_type; end
try handles.segment_mode = segment_mode; end
try handles.Seg_thres = Seg_thres; end
try handles.fish_dir = fish_dir; end
try handles.outfile_prefix_seg = outfile_prefix_seg; end
try handles.max_int_thres = max_int_thres; end
try handles.max_int_spot_th = max_int_spot_th; end
try handles.min_nucleus_size = min_nucleus_size; end
try handles.max_nucleus_size = max_nucleus_size; end
try handles.min_cell_size = min_cell_size; end
try handles.max_cell_size = max_cell_size; end
try handles.yeast_seg = yeast_seg; end
try handles.input_name = get(handles.Load_In_edit,'String'); end
try handles.Diffstack = Diffstack; catch Diffstack = 0; handles.Diffstack = Diffstack; end
try handles.TRANS_prefix = TRANS_prefix; catch TRANS_prefix = ''; handles.TRANS_prefix = TRANS_prefix; end
try handles.DAPI_prefix = DAPI_prefix; catch DAPI_prefix = ''; handles.DAPI_prefix = TRANS_prefix; end
%% Populate everything on the form with updated variables
% set(handles.Beg_thres_edit,'Value',handles.thA(1)); % Thresholds for RNA spot intensities
% set(handles.Int_edit,'Value',handles.thA(2));
% set(handles.End_Thres_edit,'Value',handles.thA(3));
%set(handles.handles.ths  % Thresholds for RNA spot intensities % YCY7 AF700 CY5 AF594 TMR YFP GFP/AF488 DAPI TRANS defines which fluorecnet dyes have been imaged 
% set(handles.DAPI_checkbox,'Value',handles.Ych(2)) ;
% set(handles.CY5_checkbox,'Value',handles.Ych(3));
% set(handles.AF594_checkbox,'Value',handles.Ych(4)) ;
% set(handles.TMR_checkbox,'Value',handles.Ych(5)) ;
% set(handles.YFP_checkbox,'Value',handles.Ych(6)) ;
% set(handles.GFP_checkbox,'Value',handles.Ych(7)) ;
% set(handles.AF700_checkbox,'Value',handles.Ych(8)); 
% set(handles.TRANS_checkbox,'Value',handles.Ych(9)) ;
set(handles.Diff_Stack_checkbox,'Value',Diffstack)
%%%Make changes depending on Diffstack (whether DAPI and TRANS images are
%%%in separate stacks)
if handles.Diffstack == 0
    set(handles.Mark_im_prefix_edit, 'enable', 'off')
    set(handles.Bound_im_prefix_edit, 'enable', 'off')
    set(handles.text18, 'enable', 'off')
    set(handles.text17, 'enable', 'off')  
    set(handles.Load_Sample_pushbutton, 'enable', 'on')
    set(handles.Boundary_slice_edit, 'enable', 'on') 
    set(handles.text19, 'enable', 'on') 
    set(handles.Marker_slice_edit, 'enable', 'on') 
    set(handles.text20, 'enable', 'on') 
    set(handles.Total_channel_edit, 'enable', 'on') 
    set(handles.text21, 'enable', 'on')       
elseif handles.Diffstack == 1
    set(handles.Mark_im_prefix_edit, 'enable', 'on')
    set(handles.Bound_im_prefix_edit, 'enable', 'on')
    set(handles.text18, 'enable', 'on')
    set(handles.text17, 'enable', 'on')
    set(handles.Load_Sample_pushbutton, 'enable', 'off')
    set(handles.Boundary_slice_edit, 'enable', 'off') 
    set(handles.text19, 'enable', 'off') 
    set(handles.Marker_slice_edit, 'enable', 'off') 
    set(handles.text20, 'enable', 'off') 
    set(handles.Total_channel_edit, 'enable', 'off') 
    set(handles.text21, 'enable', 'off')      
end
set(handles.Mark_im_prefix_edit,'String',DAPI_prefix)
set(handles.Bound_im_prefix_edit,'String',TRANS_prefix)
set(handles.Low_Slice_edit,'String',num2str(handles.img_stacks(1)));
set(handles.High_Slice_edit,'String',num2str(handles.img_stacks(2)));
set(handles.ManTh_checkbox,'Value',handles.ManTh);
%set(handles.handles.Yth
% set(handles.RNA_Dir_edit,'String',handles.outfile_prefix_RNA);
%set(handles.handles.ManTh
set(handles.Ywin_checkbox,'Value',handles.Ywin);
set(handles.Yim_checkbox,'Value',handles.Yim);
% handles.contents = cellstr(get(handles.File_Type_popupmenu,'String')); %returns File_Type_popupmenu contents as cell array
% temp_ind = 1;    %Will have the index of the file type for showing it. but defaults to "File type"
% for i = 1:size(handles.contents)
%     if strcmp(handles.file_type,handles.contents{i})    %compare file type to possible ones listed
%         temp_ind = i;
%     end
% end
% set(handles.File_Type_popupmenu,'Value',temp_ind); %set the file type in the dropdown menu to the correct one  
%set(handles.handles.segment_mode
set(handles.Single_Thres_checkbox,'Value',handles.Seg_thres);
set(handles.Im_Dir_edit,'String',handles.fish_dir);
set(handles.Seg_Dir_edit,'String',handles.outfile_prefix_seg);
% set(handles.Max_Int_Thres_checkbox,'Value',handles.max_int_thres);
% set(handles.Max_Int_Spot_Thres_checkbox,'Value',handles.max_int_spot_th);
set(handles.Min_Nuc_edit,'String',num2str(handles.min_nucleus_size));
set(handles.Max_Nuc_edit,'String',num2str(handles.max_nucleus_size));
set(handles.Min_Cell_edit,'String',num2str(handles.min_cell_size));
set(handles.Max_Cell_edit,'String',num2str(handles.max_cell_size));
set(handles.Total_channel_edit,'String',num2str(handles.Ych(1)));
set(handles.Boundary_slice_edit,'String',num2str(handles.Ych(3)));
set(handles.Marker_slice_edit,'String',num2str(handles.Ych(2)));
%handles.yeast_seg
set(handles.Input_Name_edit,'String',get(handles.Load_In_edit,'String'));
% [images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin) ;   %Gets images if not Diffstack
% handles.images1 = images1;
% handles.im_prefixes = im_prefixes;
% if handles.Diffstack
%  [images1,folders1] = Get_Dir_All(handles.fish_dir,handles.Ywin);
%  clear images2   
%  images2 = {};
%  for i = 1:size(images1,2)
%      if Ywin
%      beginning1 = strfind(images1{i},'\');
%      else
%          beginning1 = strfind(images1{i},'/');
%      end
%      beginning1 = beginning1(end)+1;    %Find the part that is not directory
%      if strfind(images1{i},'_MMStack')
%          end1 = strfind(images1{i},'_MMStack');
%      else
%          end1 = strfind(images1{i},'.');
%      end
%      end1 = end1(1)-1;
%      images2{i} = images1{i}(beginning1:end1);
%  end 
%   handles.images1 = images1;
%   handles.images2 = images2;
%   counter_DAPI = 1;
%   counter_TRANS = 1; 
%   handles.DAPI_prefix
%   for i = 1:size(images2,2)
%       if strfind(images2{i},handles.DAPI_prefix)
%           handles.images_DAPI{counter_DAPI} = images1{i};
%           counter_DAPI = counter_DAPI+1;
%       elseif strfind(images2{i},handles.TRANS_prefix)
%           handles.images_TRANS{counter_TRANS} = images1{i};
%           counter_TRANS = counter_TRANS+1;
%       end
%   end 
% else  
% [images1,im_prefixes] = Get_Dir_ImSt(handles.fish_dir,handles.Ywin);                                           %Finds the file names within the fish directory
%     handles.images1 = images1;
%     handles.im_prefixes = im_prefixes;
% end
[images1,folders1] = Get_Dir_All(handles.fish_dir,handles.Ywin);
 clear images2 im_prefixes  
 images2 = {};          %This will have the file names (not the entire directory and not the type of file)
 for i = 1:size(images1,2)
     if Ywin                     %The type of slash depends on whether it is windows
     beginning1 = strfind(images1{i},'\');      
     else
         beginning1 = strfind(images1{i},'/');
     end
     beginning1 = beginning1(end)+1;    %Find the part that is not directory
      if strfind(images1{i},'_MMStack')
         end1 = strfind(images1{i},'_MMStack'); %look for MMStack since it is often at the end of the file names
     elseif strfind(images1{i},'_MMImages')
         end1 = strfind(images1{i},'_MMImages'); %look for MMImages since it is often at the end of the file names
     elseif strfind(images1{i},'.ome')
         end1 = strfind(images1{i},'.ome');        %Look for a period otherwise as the end of the file name
     else
         end1 = strfind(images1{i},'.');        %Look for a period otherwise as the end of the file name
     end
     end1 = end1(1)-1;
     images2{i} = images1{i}(beginning1:end1);  %Populate the cell array with the file names
 end 
 images2';  %For visualizing all files extracted (remove semicolon)
  handles.images1 = images1;
  handles.images2 = images2;
if handles.Diffstack
    clear handles.images_DAPI handles.images_TRANS
    handles.images_DAPI = {};
    handles.images_DAPI = {};
  counter_DAPI = 1;             %counter for marker images found +1 (used to populate both DAPI and TRANS cell arrays)
  counter_TRANS = 0;            %counter for TRANS images found (used to determine if TRANS prefix was correct)
  handles.DAPI_prefix
  for i = 1:size(images2,2) 
      if strfind(images2{i},handles.DAPI_prefix) == 1
          handles.images_DAPI{counter_DAPI} = images1{i};
          im_prefixes{counter_DAPI} =  images2{i}(size(handles.DAPI_prefix,2)+1:end);   %populate imprefixes (used for later segmentation file names)
          for j = 1:size(images2,2)
              if strfind(images2{j},handles.TRANS_prefix) == 1 &...         %Look for TRANS prefix at beginning
                      strcmp(images2{j}(size(handles.TRANS_prefix,2)+1:end),... %Look for the rest of image name to match between TRANS
                      images2{i}(size(handles.DAPI_prefix,2)+1:end))        %and DAPI
                      handles.images_TRANS{counter_DAPI} = images1{j}; 
              end
          end
          counter_TRANS = counter_TRANS+1;                                  %Add to counter (even if match not found)
          counter_DAPI = counter_DAPI+1;                                    %Add to counter (even if match not found)
       if strfind(images2{i},handles.TRANS_prefix) %searching for TRANS prefix
%            handles.images_TRANS{counter_TRANS} = images1{i};
           counter_TRANS = counter_TRANS+1;
       end
      end
  end
  if counter_TRANS == 0
      msgbox('There were no images found with the boundary prefix')
  end
  if isempty(handles.images_DAPI)
      msgbox('There were no images found with the marker prefix')     
  elseif isempty(handles.images_TRANS)
      msgbox('There were no images paired with the prefixes specified')
  end
else
%[images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                %Finds the file names within the fish directory (old way)
im_prefixes = images2;
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in Load_Sample_pushbutton.
function Load_Sample_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Sample_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[file,path] = uigetfile;
f = msgbox('Loading Data');
imnum1 = 1;
fail1 = 1;
handles.images1{1}
[stack, img_read] = tiffread2(handles.images1{1});
handles.sample_file = handles.images1{1};
 Ych = handles.Ych
%     chi = find(handles.Ych == 1);                                                       % Determines how many channels have been imaged
%     ch = size(chi,2);
%in Ych, first is the total number of channels, second is the marker slice, and third is the boundary slice
ch = Ych(1)          
ImSt = img_read/ch

        if Ych(2)                                                      % load DAPI    CHANGED TO 7 TEMPORARILY BECAUSE OF SOMETHING STRANGE                             
            ij = Ych(2);                                                  %CHANGED TO 7 TEMPORARILY BECAUSE OF SOMETHING STRANGE
            DAPIim = [ij:ch:img_read];
            im_size = size(stack(1,DAPIim(1)).data);                        %BK 5/15/15 
            DAPI_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/15/15 resets the DAPI image so previous ones do not get incorporated
            for i = 1:ImSt;
                DAPI_ims(:,:,i) = stack(1,DAPIim(i)).data;
            end;
            handles.DAPI_ims = DAPI_ims;
            handles.nuc_border = zeros(size(DAPI_ims));
            handles.cell_border = zeros(size(DAPI_ims,1),size(DAPI_ims,2));
        else
            im_size = size(stack(1,1).data);
            DAPI_ims = zeros(im_size(1),im_size(2),ImSt);
            handles.DAPI_ims = DAPI_ims;
            handles.nuc_border = zeros(size(DAPI_ims));
            handles.cell_border = zeros(size(DAPI_ims,1),size(DAPI_ims,2));
        end;
        if Ych(3)                                                     % load TRANS
            ij = Ych(3);
            TRANSim = [ij:ch:img_read];
            im_size = size(stack(1,TRANSim(1)).data);                        %BK 5/15/15
            TRANS_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 5/15/15 resets the TRANS image so previous ones do not get incorporated
            for i = 1:ImSt
                TRANS_ims(:,:,i) = stack(1,TRANSim(i)).data;
            end
        else
            im_size = size(stack(1,1).data);
            TRANS_ims = zeros(im_size(1),im_size(2),ImSt);
            handles.TRANS_ims = TRANS_ims;
            handles.nuc_border = zeros(size(TRANS_ims));
            handles.cell_border = zeros(size(TRANS_ims,1),size(TRANS_ims,2));
        end;
        handles.TRANS_ims = TRANS_ims;
        set(handles.Slice_viewed_boundary_slider,'Min',1)
        set(handles.Slice_viewed_boundary_slider,'Max',ImSt)
        set(handles.Slice_viewed_DAPI_slider,'Min',1)
        set(handles.Slice_viewed_DAPI_slider,'Max',ImSt)
        set(handles.Slice_viewed_DAPI_slider,'SliderStep', [1/ImSt, 0.1])
        set(handles.Slice_viewed_boundary_slider,'SliderStep', [1/ImSt, 0.1])
        set(handles.Slice_viewed_DAPI_slider,'Value',round(ImSt/2))
        set(handles.Slice_viewed_boundary_slider,'Value',round(ImSt/2))
        set(handles.DAPI_slice_text,'String',['Slice Number ' num2str(round(ImSt/2))])
        set(handles.TRANS_slice_text,'String',['Slice Number ' num2str(round(ImSt/2))])
        handles.DAPI_slice = round(ImSt/2);
        handles.TRANS_slice = round(ImSt/2);
        Mid_DAPI = max(DAPI_ims(:,:,round(ImSt/2)),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[0 1]);
            Trans1 = max(TRANS_ims(:,:,round(ImSt/2)),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[0 1]);
            R = Trans3;
            B = Nuc1+Trans3;
            G = Trans3;
            RGB = cat(3,R,G,B);     
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])            
            guidata(hObject, handles);
            set(handles.Segment_sample_pushbutton,'Enable','on')
            set(handles.Segment_fourth_pushbutton,'Enable','on')
fail1 = 0;
            %removeToolbarExplorationButtons(gcf)
            %set(handles.Sample_Image_axes,'toolbar','figure');




function Beg_thres_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Beg_thres_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.thA(1) = str2num(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Beg_thres_edit as text
%        str2double(get(hObject,'String')) returns contents of Beg_thres_edit as a double


% --- Executes during object creation, after setting all properties.
function Beg_thres_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beg_thres_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Int_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Int_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.thA(2) = str2num(get(hObject,'String')); % Thresholds for RNA spot intensities
% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Int_edit as text
%        str2double(get(hObject,'String')) returns contents of Int_edit as a double


% --- Executes during object creation, after setting all properties.
function Int_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Int_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function End_Thres_edit_Callback(hObject, eventdata, handles)
% hObject    handle to End_Thres_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.thA(3) = str2num(get(hObject,'String')); % Thresholds for RNA spot intensities
% Update handles structure
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of End_Thres_edit as text
%        str2double(get(hObject,'String')) returns contents of End_Thres_edit as a double


% --- Executes during object creation, after setting all properties.
function End_Thres_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to End_Thres_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RNA_Dir_pushbutton.
function RNA_Dir_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to RNA_Dir_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.outfile_prefix_RNA = uigetdir();
% Update handles structure
guidata(hObject, handles);
set(handles.RNA_Dir_edit,'String',handles.outfile_prefix_RNA);



function RNA_Dir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to RNA_Dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RNA_Dir_edit as text
%        str2double(get(hObject,'String')) returns contents of RNA_Dir_edit as a double
handles.outfile_prefix_RNA = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function RNA_Dir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RNA_Dir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Max_Int_Thres_checkbox.
function Max_Int_Thres_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Max_Int_Thres_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.max_int_thres = get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Max_Int_Thres_checkbox


% --- Executes on button press in Draw_Ellipse_pushbutton.
function Draw_Ellipse_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Draw_Ellipse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.Sample_Image_axes);
% imshow(RGB,[]);
% ax = gca;
% ax.Toolbar.Visible = 'on';
h = imellipse;
%// Create a logical mask
roi = createMask(h);

%// Sum the values equal to 1;
Area = sum(roi(:));
set(handles.Ellipse_Area_edit,'String',num2str(Area));
handles.area = Area;
guidata(hObject, handles);



function Ellipse_Area_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Ellipse_Area_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ellipse_Area_edit as text
%        str2double(get(hObject,'String')) returns contents of Ellipse_Area_edit as a double


% --- Executes during object creation, after setting all properties.
function Ellipse_Area_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ellipse_Area_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Select_Sample_pushbutton.
function Select_Sample_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Sample_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
images1 = handles.images1;
if handles.Diffstack
    [indx1,tf1] = listdlg('ListString',images1,'ListSize',[1000 500],'SelectionMode','single',...
        'PromptString','Select the marker (DAPI) image to process') % creates a modal dialog box that allows the user to select one or more items from the specified list.   
     [indx2,tf2] = listdlg('ListString',images1,'ListSize',[1000 500],'SelectionMode','single',...
        'PromptString','Select the boundary(TRANS) image to process') % creates a modal dialog box that allows the user to select one or more items from the specified list.   
    if tf1 & tf2
        handles.indx1 = indx1;
        handles.indx2 = indx2;
        f = msgbox('Loading Data');
        [stack, img_read] = tiffread2([images1{indx1}]);
        im_size = size(stack(1,1).data);                        %BK 5/15/15
        DAPI_ims = NaN(im_size(1),im_size(2),img_read);                     %BK 5/15/15 resets the DAPI image so previous ones do not get incorporated
        for i = 1:img_read
            DAPI_ims(:,:,i) = stack(1,i).data;
        end
        clear stack
        [stack, img_read] = tiffread2([images1{indx2}]);
        im_size = size(stack(1,1).data);                        %BK 5/15/15
        TRANS_ims = NaN(im_size(1),im_size(2),img_read);                     %BK 5/15/15 resets the DAPI image so previous ones do not get incorporated
        for i = 1:img_read
            TRANS_ims(:,:,i) = stack(1,i).data;
        end
        ImSt = img_read;
        clear stack
    end
else
    [file,path] = uigetfile;
    handles.sample_file = [path file];
    f = msgbox('Loading Data');
    [stack, img_read] = tiffread2([path file]);
    Ych = handles.Ych;
    ch = Ych(1);
    ImSt = img_read/ch;
    if Ych(2)                                                      % load DAPI    CHANGED TO 7 TEMPORARILY BECAUSE OF SOMETHING STRANGE
        ij = Ych(2);                                                  %CHANGED TO 7 TEMPORARILY BECAUSE OF SOMETHING STRANGE
        DAPIim = [ij:ch:img_read];
        im_size = size(stack(1,DAPIim(1)).data);                        %BK 5/15/15
        DAPI_ims = NaN(im_size(1),im_size(2),ImSt);                     %BK 5/15/15 resets the DAPI image so previous ones do not get incorporated
        for i = 1:ImSt;
            DAPI_ims(:,:,i) = stack(1,DAPIim(i)).data;
        end;
        handles.nuc_border = zeros(size(DAPI_ims));
        handles.cell_border = zeros(size(DAPI_ims,1),size(DAPI_ims,2));
    else
    end;
    
    if Ych(3)                                                     % load TRANS
        ij = Ych(3);
        TRANSim = [ij:ch:img_read];
        im_size = size(stack(1,TRANSim(1)).data);                        %BK 5/15/15
        TRANS_ims = NaN(im_size(1),im_size(2),ImSt);                    %BK 5/15/15 resets the TRANS image so previous ones do not get incorporated
        for i = 1:ImSt;
            TRANS_ims(:,:,i) = stack(1,TRANSim(i)).data;
        end;
    else
    end;
end
handles.DAPI_ims = DAPI_ims;
handles.TRANS_ims = TRANS_ims;
set(handles.Slice_viewed_boundary_slider,'Min',1)
set(handles.Slice_viewed_boundary_slider,'Max',ImSt)
set(handles.Slice_viewed_DAPI_slider,'Min',1)
set(handles.Slice_viewed_DAPI_slider,'Max',ImSt)
set(handles.Slice_viewed_DAPI_slider,'SliderStep', [1/ImSt, 0.1])
set(handles.Slice_viewed_boundary_slider,'SliderStep', [1/ImSt, 0.1])
set(handles.Slice_viewed_DAPI_slider,'Value',round(ImSt/2))
set(handles.Slice_viewed_boundary_slider,'Value',round(ImSt/2))
handles.DAPI_slice = round(ImSt/2);
handles.TRANS_slice = round(ImSt/2);
Mid_DAPI = max(DAPI_ims(:,:,round(ImSt/2)),[],3);
Nuc1 = Mid_DAPI/max(DAPI_ims(:));
Nuc2 = imadjust(Nuc1,[0 1]);
Trans1 = max(TRANS_ims(:,:,round(ImSt/2)),[],3);
Trans2 = Trans1/max(TRANS_ims(:));
Trans3 = imadjust(Trans2,[0 1]);
R = Trans3;
B = Nuc1+Trans3;
G = Trans3;
RGB = cat(3,R,G,B);
axes(handles.Sample_Image_axes);
imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
ax = gca;
ax.Toolbar.Visible = 'on';
set(handles.Sample_Image_axes,'xticklabel',[])
set(handles.Sample_Image_axes,'xtick',[])
set(handles.Sample_Image_axes,'yticklabel',[])
set(handles.Sample_Image_axes,'ytick',[])
set(handles.Segment_sample_pushbutton,'Enable','on')
set(handles.Segment_fourth_pushbutton,'Enable','on')
guidata(hObject, handles);
%removeToolbarExplorationButtons(gcf)
%set(handles.Sample_Image_axes,'toolbar','figure');


% --- Executes on button press in ManTh_checkbox.
function ManTh_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to ManTh_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ManTh = get(hObject,'Value');
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of ManTh_checkbox



function Bound_im_prefix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Bound_im_prefix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.TRANS_prefix = get(hObject,'String');
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Bound_im_prefix_edit as text
%        str2double(get(hObject,'String')) returns contents of Bound_im_prefix_edit as a double


% --- Executes during object creation, after setting all properties.
function Bound_im_prefix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bound_im_prefix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Mark_im_prefix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Mark_im_prefix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.DAPI_prefix = get(hObject,'String');
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of Mark_im_prefix_edit as text
%        str2double(get(hObject,'String')) returns contents of Mark_im_prefix_edit as a double


% --- Executes during object creation, after setting all properties.
function Mark_im_prefix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Mark_im_prefix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Diff_Stack_checkbox.
function Diff_Stack_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Diff_Stack_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Diffstack = get(hObject,'Value');
if handles.Diffstack == 0
    set(handles.Mark_im_prefix_edit, 'enable', 'off')
    set(handles.Bound_im_prefix_edit, 'enable', 'off')
    set(handles.text18, 'enable', 'off')
    set(handles.text17, 'enable', 'off')  
    set(handles.Load_Sample_pushbutton, 'enable', 'on')
    set(handles.Boundary_slice_edit, 'enable', 'on') 
    set(handles.text19, 'enable', 'on') 
    set(handles.Marker_slice_edit, 'enable', 'on') 
    set(handles.text20, 'enable', 'on') 
    set(handles.Total_channel_edit, 'enable', 'on') 
    set(handles.text21, 'enable', 'on')       
elseif handles.Diffstack == 1
    set(handles.Mark_im_prefix_edit, 'enable', 'on')
    set(handles.Bound_im_prefix_edit, 'enable', 'on')
    set(handles.text18, 'enable', 'on')
    set(handles.text17, 'enable', 'on')
    set(handles.Load_Sample_pushbutton, 'enable', 'off')
    set(handles.Boundary_slice_edit, 'enable', 'off') 
    set(handles.text19, 'enable', 'off') 
    set(handles.Marker_slice_edit, 'enable', 'off') 
    set(handles.text20, 'enable', 'off') 
    set(handles.Total_channel_edit, 'enable', 'off') 
    set(handles.text21, 'enable', 'off')      
end
% Update handles structure
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Diff_Stack_checkbox



function Boundary_slice_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Boundary_slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%first is the total number of channels, second is the marker slice, and third is the boundary slice
Ych = handles.Ych;
Ych(3) = str2num(get(hObject,'String'));
handles.Ych = Ych;
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of Boundary_slice_edit as text
%        str2double(get(hObject,'String')) returns contents of Boundary_slice_edit as a double


% --- Executes during object creation, after setting all properties.
function Boundary_slice_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Boundary_slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Marker_slice_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Marker_slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%first is the total number of channels, second is the marker slice, and third is the boundary slice
Ych = handles.Ych;
Ych(2) = str2num(get(hObject,'String'));
handles.Ych = Ych;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Marker_slice_edit as text
%        str2double(get(hObject,'String')) returns contents of Marker_slice_edit as a double


% --- Executes during object creation, after setting all properties.
function Marker_slice_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Marker_slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Total_channel_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Total_channel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%first is the total number of channels, second is the marker slice, and third is the boundary slice
Ych = handles.Ych;
Ych(1) = str2num(get(hObject,'String'));
handles.Ych = Ych;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Total_channel_edit as text
%        str2double(get(hObject,'String')) returns contents of Total_channel_edit as a double


% --- Executes during object creation, after setting all properties.
function Total_channel_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Total_channel_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function Min_brightness_DAPI_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Min_brightness_DAPI_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
val= get(hObject,'Value');
%'This happened'
handles.DAPI_min = val;
if handles.DAPI_max <= val
    handles.DAPI_max = val+.01
    set(handles.Max_brightness_DAPI_slider,'Value',handles.DAPI_max)
end
%hObject.Value=val;
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                  CellBorder1 = imclose(CellBorder1,se); 
                  if false %handles.thicken_num >1
                     CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); 
                  end
            else
                CellBorder1 = zeros(size(Trans3));
            end
            R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])
            guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Min_brightness_DAPI_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_brightness_DAPI_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Max_brightness_DAPI_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Max_brightness_DAPI_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
val= get(hObject,'Value');
handles.DAPI_max = val;
if handles.DAPI_min >= val;
    handles.DAPI_min = val-.01
    set(handles.Min_brightness_DAPI_slider,'Value',handles.DAPI_min)
end
%hObject.Value=val;
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                  CellBorder1 = imclose(CellBorder1,se);
                  if false %handles.thicken_num >1
                      CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se);
                  end
            else
                CellBorder1 = zeros(size(Trans3));
            end
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);     
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)  
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])            
            guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Max_brightness_DAPI_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_brightness_DAPI_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Slice_viewed_DAPI_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Slice_viewed_DAPI_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
val= round(get(hObject,'Value'))
handles.DAPI_slice = val;
hObject.Value=val;
set(handles.DAPI_slice_text,'String',['Slice Number ' num2str(val)])
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                  CellBorder1 = imclose(CellBorder1,se); 
                  if false %handles.thicken_num >1
                      CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se);
                  end
            else
                CellBorder1 = zeros(size(Trans3));
            end
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);     
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])            
            guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Slice_viewed_DAPI_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slice_viewed_DAPI_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Min_brightness_boundary_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Min_brightness_boundary_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
val= get(hObject,'Value');
handles.TRANS_min = val;
if handles.TRANS_max <= val;
    handles.TRANS_max = val+.01
    set(handles.Max_brightness_boundary_slider,'Value',handles.TRANS_max)
end
%hObject.Value=val;
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                  CellBorder1 = imclose(CellBorder1,se); 
                  if false %handles.thicken_num >1
                  CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); 
                  end
            else
                CellBorder1 = zeros(size(Trans3));
            end
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);     
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])            
            guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Min_brightness_boundary_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Min_brightness_boundary_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Max_brightness_boundary_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Max_brightness_boundary_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
val= get(hObject,'Value');
handles.TRANS_max = val;
if handles.TRANS_min >= val;
    handles.TRANS_min = val-.01
    set(handles.Min_brightness_TRANS_slider,'Value',handles.TRANS_min)
end
%hObject.Value=val;
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                  CellBorder1 = imclose(CellBorder1,se); 
                  if false %handles.thicken_num >1
                   CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); 
                  end
                  else
                CellBorder1 = zeros(size(Trans3));
            end
            R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);     
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])            
            guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Max_brightness_boundary_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Max_brightness_boundary_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function Slice_viewed_boundary_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Slice_viewed_boundary_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
val= round(get(hObject,'Value'));
handles.TRANS_slice = val;
hObject.Value=val;
set(handles.TRANS_slice_text,'String',['Slice Number ' num2str(val)])
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                  CellBorder1 = imclose(CellBorder1,se); 
                  if false %handles.thicken_num >1
                  CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); 
                  end
             else
                CellBorder1 = zeros(size(Trans3));
            end
            R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);     
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')            
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])            
            guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Slice_viewed_boundary_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slice_viewed_boundary_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on mouse press over axes background.
function Sample_Image_axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Sample_Image_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function Sample_Image_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sample_Image_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: place code in OpeningFcn to populate Sample_Image_axes


% --- Executes on slider movement.
function Bound_thick_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Bound_thick_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
val= round(get(hObject,'Value'));
handles.thicken_num = val;
hObject.Value=val;
%set(handles.TRANS_slice_text,'String',['Slice Number ' num2str(val)])
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);                
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                  CellBorder1 = imclose(CellBorder1,se); 
                  if false %handles.thicken_num >1
                  CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se);             
                  end
            else
                  CellBorder1 = zeros(size(Trans3));
            end
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);     
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])            
            guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function Bound_thick_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bound_thick_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Segment_sample_pushbutton.
function Segment_sample_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Segment_sample_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;
ths = handles.ths;
Yth = handles.Yth;
Ych = handles.Ych;
ManTh = handles.ManTh;
Ywin = handles.Ywin;
Yim = handles.Yim;
file_type = handles.file_type;
segment_mode = handles.segment_mode;
Seg_thres = handles.Seg_thres;
fish_dir = handles.fish_dir;
outfile_prefix_seg = handles.outfile_prefix_seg;
node_num = handles.node_num;
max_int_thres = handles.max_int_thres;
max_int_spot_th = handles.max_int_spot_th;
min_nucleus_size = handles.min_nucleus_size;
max_nucleus_size = handles.max_nucleus_size;
min_cell_size = handles.min_cell_size;
max_cell_size = handles.max_cell_size;
Diffstack=handles.Diffstack;
img_stacks = handles.img_stacks;
yeast_seg = handles.yeast_seg;
if Diffstack
    images1 = handles.images1;   
    DAPI_images = {images1{handles.indx1}};
    TRANS_images = {images1{handles.indx2}};
    images1 = DAPI_images;
else
    images1 = {handles.sample_file};
    images1
    DAPI_images = {};
    TRANS_images = {};
end
im_prefixes = {'Sample_Seg'};
node_num = handles.node_num;
% [images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                                     %Finds the file names within the fish directory
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
dxy1 = round(sqrt(max_nucleus_size/pi)); 
A1_segment_predefined_variables_streamlined_generalized(Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_stacks,yeast_seg,images1,im_prefixes,dxy1,Diffstack,DAPI_images,TRANS_images)
%%% Display segmented boundaries
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
load([outfile_prefix_seg 'Lab_Sample_Seg'],'cells')
load([outfile_prefix_seg 'nuclei_Sample_Seg'],'Label_mid')
handles.nuc_border = Label_mid;
handles.cell_border = cells;
NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
se = strel('disk',6);
NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
se = strel('disk',6);
  CellBorder1 = imclose(CellBorder1,se); 
  if false %handles.thicken_num >1
  CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); 
  end
  DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
Nuc1 = Mid_DAPI/max(DAPI_ims(:));
Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
Trans2 = Trans1/max(TRANS_ims(:));
Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
'size Trans3'
size(Trans3)
'size NucBorder1'
size(NucBorder1)
'size Nuc2'
size(Nuc2)
'size CellBorder1'
size(CellBorder1)
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
RGB = cat(3,R,G,B);
axes(handles.Sample_Image_axes);
imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
ax = gca;
ax.Toolbar.Visible = 'on';
set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
set(handles.Sample_Image_axes,'xticklabel',[])
set(handles.Sample_Image_axes,'xtick',[])
set(handles.Sample_Image_axes,'yticklabel',[])
set(handles.Sample_Image_axes,'ytick',[])
set(handles.Bound_thick_slider,'SliderStep', [1/11, 0.1])
set(handles.Bound_thick_slider,'Min', 0)
set(handles.Bound_thick_slider,'Max', 10)
set(handles.Bound_thick_slider,'Enable', 'On')
set(handles.Bound_thick_text,'Enable', 'On')
set(handles.Nuc_border_checkbox,'Enable', 'On')
set(handles.Cell_border_checkbox,'Enable', 'On')
set(handles.Nuc_border_checkbox,'Value', 1)
set(handles.Cell_border_checkbox,'Value', 1)
guidata(hObject, handles);


% --- Executes on button press in Nuc_border_checkbox.
function Nuc_border_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Nuc_border_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.show_nuc_bord = get(hObject,'Value');
guidata(hObject, handles);
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
%hObject.Value=val;
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                  CellBorder1 = imclose(CellBorder1,se);
                  if false %handles.thicken_num >1
                      CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); 
                  end
            else
                CellBorder1 = zeros(size(Trans3));
            end
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])
            guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Nuc_border_checkbox


% --- Executes on button press in Cell_border_checkbox.
function Cell_border_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to Cell_border_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.show_cell_bord = get(hObject,'Value');
guidata(hObject, handles);
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
%hObject.Value=val;
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                 CellBorder1 = imclose(CellBorder1,se); 
                 if false %handles.thicken_num >1
                     CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se);
                 end
            else          
                CellBorder1 = zeros(size(Trans3));
            end
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);
            axes(handles.Sample_Image_axes);
            imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
            ax = gca;
            ax.Toolbar.Visible = 'on';
            set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
            set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
            set(handles.Sample_Image_axes,'xticklabel',[])
            set(handles.Sample_Image_axes,'xtick',[])
            set(handles.Sample_Image_axes,'yticklabel',[])
            set(handles.Sample_Image_axes,'ytick',[])
            guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of Cell_border_checkbox


% --- Executes on button press in Save_image_pushbutton.
function Save_image_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Save_image_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% imsave(handles.Sample_Image_axes);
x_vals = get(handles.Sample_Image_axes,'Xlim')
y_vals = get(handles.Sample_Image_axes,'Ylim')
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
x_vals(x_vals>size(DAPI_ims,1)) = size(DAPI_ims,1);
y_vals(y_vals>size(DAPI_ims,1)) = size(DAPI_ims,1);
x_vals(x_vals<1) = 1;
y_vals(y_vals<1) = 1;
% for z = 1:size(DAPI_ims,3)%[1,7,14,22]
    z = handles.DAPI_slice; 
    Mid_DAPI = max(DAPI_ims(y_vals(1):y_vals(2),x_vals(1):x_vals(2),z),[],3);
    Nuc1 = Mid_DAPI/max(DAPI_ims(:));
    Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
    z = handles.TRANS_slice; 
    Trans1 = max(TRANS_ims(y_vals(1):y_vals(2),x_vals(1):x_vals(2),z),[],3);
    Trans2 = Trans1/max(TRANS_ims(:));
    Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
    if handles.show_nuc_bord == 1;
        z = handles.DAPI_slice; 
        NucBorder1 = bwmorph(handles.nuc_border(y_vals(1):y_vals(2),x_vals(1):x_vals(2),z),'remove')*100000;
        NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
       se = strel('disk',4);
        NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
    else
        NucBorder1 = zeros(size(Trans3));
    end
    if handles.show_cell_bord == 1;
        z = handles.TRANS_slice; 
        CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
        CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
       se = strel('disk',4);
         CellBorder1 = imclose(CellBorder1,se); 
         if false %handles.thicken_num >1
             CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se);
         end
         else
        CellBorder1 = zeros(size(Trans3));
    end
    R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
    RGB = cat(3,R,G,B);
%     imwrite(RGB, ['Sample_ImageStackSlice' num2str(z) '.tiff'],'Compression','none')
if not(isempty(handles.Im_name))
imwrite(RGB, [handles.Im_name '.tiff'],'Compression','none')
else
    imwrite(RGB, ['Saved_Image.tiff'],'Compression','none')
end
% end
%             Tiff(RGB)
%             save('Sample_RGB','RGB','-v7.3')
%             imwrite('Sample_Image_stack','RGB','.tiff')
%             figure();
%             imshow(RGB,[]);
% F = getframe(handles.Sample_Image_axes);
% saveas(handles.Sample_Image_axes, [handles.Im_name '.eps'])
% saveas(handles.Sample_Image_axes,[handles.Im_name '.fig'])
% Image = frame2im(F);
% imwrite(Image, [handles.Im_name '.jpg'])
% imwrite(Image, [handles.Im_name '.png'])
% imwrite(Image, [handles.Im_name '.tiff'])
% imwrite(Image, 'Image.eps')


% --- Executes on selection change in Presets_popupmenu.
function Presets_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to Presets_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
preset1 = contents{get(hObject,'Value')};
if strcmp(preset1,'100x S. Cerevisiae')
    set(handles.Min_Nuc_edit,'String','350');
    set(handles.Max_Nuc_edit,'String','1500');
    set(handles.Min_Cell_edit,'String','2000');
    set(handles.Max_Cell_edit,'String','15000');
    set(handles.Low_Slice_edit,'String','15');
    set(handles.High_Slice_edit,'String','26');
elseif strcmp(preset1,'100x S. Pombe')
    set(handles.Min_Nuc_edit,'String','250');
    set(handles.Max_Nuc_edit,'String','800');
    set(handles.Min_Cell_edit,'String','3500');
    set(handles.Max_Cell_edit,'String','20000');
    set(handles.Low_Slice_edit,'String','15');
    set(handles.High_Slice_edit,'String','19');
elseif strcmp(preset1,'20x S. Pombe')
    set(handles.Min_Nuc_edit,'String','20');
    set(handles.Max_Nuc_edit,'String','35');
    set(handles.Min_Cell_edit,'String','195');
    set(handles.Max_Cell_edit,'String','600');
    set(handles.Low_Slice_edit,'String','15');
    set(handles.High_Slice_edit,'String','17');
elseif strcmp(preset1,'100x mESCs')
    set(handles.Min_Nuc_edit,'String','10000');
    set(handles.Max_Nuc_edit,'String','35000');
    set(handles.Min_Cell_edit,'String','20000');
    set(handles.Max_Cell_edit,'String','60000');
    set(handles.Low_Slice_edit,'String','36');
    set(handles.High_Slice_edit,'String','65');
elseif strcmp(preset1,'20x mESCs')
    set(handles.Min_Nuc_edit,'String','200');
    set(handles.Max_Nuc_edit,'String','2300');
    set(handles.Min_Cell_edit,'String','300');
    set(handles.Max_Cell_edit,'String','3000');
    set(handles.Low_Slice_edit,'String','25');
    set(handles.High_Slice_edit,'String','33');
elseif strcmp(preset1,'100x Human (Jurkat)')
    set(handles.Min_Nuc_edit,'String','7000');
    set(handles.Max_Nuc_edit,'String','30000');
    set(handles.Min_Cell_edit,'String','15000');
    set(handles.Max_Cell_edit,'String','60000');
    set(handles.Low_Slice_edit,'String','21');
    set(handles.High_Slice_edit,'String','30');
elseif strcmp(preset1,'20x Human (Jurkat)')
    set(handles.Min_Nuc_edit,'String','200');
    set(handles.Max_Nuc_edit,'String','2000');
    set(handles.Min_Cell_edit,'String','300');
    set(handles.Max_Cell_edit,'String','3000');
    set(handles.Low_Slice_edit,'String','25');
    set(handles.High_Slice_edit,'String','37');
end
handles.min_nucleus_size = str2num(get(handles.Min_Nuc_edit,'String'));
handles.max_nucleus_size = str2num(get(handles.Max_Nuc_edit,'String'));
handles.min_cell_size = str2num(get(handles.Min_Cell_edit,'String'));
handles.max_cell_size = str2num(get(handles.Max_Cell_edit,'String'));
handles.img_stacks(1) =  str2num(get(handles.Low_Slice_edit,'String'));
handles.img_stacks(2) =  str2num(get(handles.High_Slice_edit,'String'));
guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns Presets_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Presets_popupmenu


% --- Executes during object creation, after setting all properties.
function Presets_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Presets_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FontSize_edit_Callback(hObject, eventdata, handles)
% hObject    handle to FontSize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_font = str2num(get(hObject,'String'));
set(handles.pushbutton19,'Fontsize',new_font);
set(handles.text50,'Fontsize',new_font);
set(handles.text51,'Fontsize',new_font);
set(handles.text49,'Fontsize',new_font);
set(handles.RNA_channel1_edit,'Fontsize',new_font);
set(handles.RNA_channel2_edit,'Fontsize',new_font);
set(handles.RNA_channel3_edit,'Fontsize',new_font);
set(handles.checkbox20,'FontSize',new_font);
set(handles.Beginning_thresh_edit,'FontSize',new_font);
set(handles.Interval_thresh_edit,'FontSize',new_font);
set(handles.End_thresh_edit,'FontSize',new_font);
set(handles.text48,'FontSize',new_font);
set(handles.Fig_pushbutton,'FontSize',new_font);
set(handles.Start_RNA_Thresholding_pushbutton,'FontSize',new_font);
set(handles.Segment_fourth_pushbutton,'FontSize',new_font);
set(handles.text46,'FontSize',new_font);
set(handles.text21,'FontSize',new_font);
set(handles.Total_channel_edit,'FontSize',new_font);
set(handles.Im_Dir_Push,'FontSize',new_font);
set(handles.Im_Dir_edit,'FontSize',new_font);
set(handles.Seg_Dir_pushbutton,'FontSize',new_font);
set(handles.Seg_Dir_edit,'FontSize',new_font);
set(handles.text9,'FontSize',new_font);
set(handles.text19,'FontSize',new_font);
set(handles.text20,'FontSize',new_font);
set(handles.Boundary_slice_edit,'FontSize',new_font);
set(handles.Marker_slice_edit,'FontSize',new_font);
set(handles.Seg_Start_pushbutton,'FontSize',new_font);
set(handles.Yim_checkbox,'FontSize',new_font);
set(handles.Ywin_checkbox,'FontSize',new_font);
set(handles.Single_Thres_checkbox,'FontSize',new_font);
% set(handles.File_Type_popupmenu,'FontSize',new_font);
set(handles.ManTh_checkbox,'FontSize',new_font);
set(handles.text17,'FontSize',new_font);
set(handles.text18,'FontSize',new_font);
set(handles.Bound_im_prefix_edit,'FontSize',new_font);
set(handles.Mark_im_prefix_edit,'FontSize',new_font);
set(handles.Diff_Stack_checkbox,'FontSize',new_font);
set(handles.Presets_popupmenu,'FontSize',new_font);
set(handles.text3,'FontSize',new_font);
set(handles.text4,'FontSize',new_font);
set(handles.text5,'FontSize',new_font);
set(handles.text6,'FontSize',new_font);
set(handles.Min_Nuc_edit,'FontSize',new_font);
set(handles.Max_Nuc_edit,'FontSize',new_font);
set(handles.Min_Cell_edit,'FontSize',new_font);
set(handles.Max_Cell_edit,'FontSize',new_font);
set(handles.text8,'FontSize',new_font);
set(handles.Low_Slice_edit,'FontSize',new_font);
set(handles.text10,'FontSize',new_font);
set(handles.High_Slice_edit,'FontSize',new_font);
set(handles.text11,'FontSize',new_font);
set(handles.Input_Name_edit,'FontSize',new_font);
set(handles.Save_In_pushbutton,'FontSize',new_font);
set(handles.Load_In_edit,'FontSize',new_font);
set(handles.Load_In_pushbutton,'FontSize',new_font);
set(handles.Load_In_Browse_pushbutton,'FontSize',new_font);
set(handles.text12,'FontSize',new_font);
set(handles.uipanel5,'FontSize',new_font);
set(handles.uipanel6,'FontSize',new_font);
set(handles.uipanel1,'FontSize',new_font);
set(handles.uipanel4,'FontSize',new_font);
set(handles.uipanel3,'FontSize',new_font);
set(handles.Load_Sample_pushbutton,'FontSize',new_font);
set(handles.Draw_Ellipse_pushbutton,'FontSize',new_font);
set(handles.Ellipse_Area_edit,'FontSize',new_font);
set(handles.text16,'FontSize',new_font);
set(handles.Select_Sample_pushbutton,'FontSize',new_font);
set(handles.Slice_viewed_DAPI_slider,'FontSize',new_font);
set(handles.Min_brightness_DAPI_slider,'FontSize',new_font);
set(handles.Max_brightness_DAPI_slider,'FontSize',new_font);
set(handles.Min_brightness_boundary_slider,'FontSize',new_font);
set(handles.Max_brightness_boundary_slider,'FontSize',new_font);
set(handles.text23,'FontSize',new_font);
set(handles.text24,'FontSize',new_font);
set(handles.text27,'FontSize',new_font);
set(handles.text28,'FontSize',new_font);
set(handles.text30,'FontSize',new_font);
set(handles.text26,'FontSize',new_font);
set(handles.Slice_viewed_boundary_slider,'FontSize',new_font);
set(handles.DAPI_slice_text,'FontSize',new_font);
set(handles.TRANS_slice_text,'FontSize',new_font);
set(handles.Bound_thick_slider,'FontSize',new_font);
set(handles.Bound_thick_text,'FontSize',new_font);
set(handles.Segment_sample_pushbutton,'FontSize',new_font);
set(handles.Nuc_border_checkbox,'FontSize',new_font);
set(handles.Cell_border_checkbox,'FontSize',new_font);
set(handles.Save_image_pushbutton,'FontSize',new_font);
set(handles.edit28,'FontSize',new_font);

guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of FontSize_edit as text
%        str2double(get(hObject,'String')) returns contents of FontSize_edit as a double


% --- Executes during object creation, after setting all properties.
function FontSize_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FontSize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit28_Callback(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Im_name = get(hObject,'String');
guidata(hObject, handles)

% Hints: get(hObject,'String') returns contents of edit28 as text
%        str2double(get(hObject,'String')) returns contents of edit28 as a double


% --- Executes during object creation, after setting all properties.
function edit28_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Fig_pushbutton.
function Fig_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Fig_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
         Nuc1 = Mid_DAPI/max(DAPI_ims(:));
            Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
            Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
            Trans2 = Trans1/max(TRANS_ims(:));
            Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
            if handles.show_nuc_bord == 1;
                NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
                NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
            else
                NucBorder1 = zeros(size(Trans3));
            end
            if handles.show_cell_bord == 1;
                CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
                CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
               se = strel('disk',4);
                 CellBorder1 = imclose(CellBorder1,se); 
                 if false %handles.thicken_num >1
                     CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se);
                 end
            else
                CellBorder1 = zeros(size(Trans3));
            end
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
            RGB = cat(3,R,G,B);  
            figure();
            imshow(RGB,[]);


% --- Executes on button press in Segment_fourth_pushbutton.
function Segment_fourth_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Segment_fourth_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;
ths = handles.ths;
Yth = handles.Yth;
Ych = handles.Ych;
ManTh = handles.ManTh;
Ywin = handles.Ywin;
Yim = handles.Yim;
file_type = handles.file_type;
segment_mode = handles.segment_mode;
Seg_thres = handles.Seg_thres;
fish_dir = handles.fish_dir;
outfile_prefix_seg = handles.outfile_prefix_seg;
node_num = handles.node_num;
max_int_thres = handles.max_int_thres;
max_int_spot_th = handles.max_int_spot_th;
min_nucleus_size = handles.min_nucleus_size;
max_nucleus_size = handles.max_nucleus_size;
min_cell_size = handles.min_cell_size;
max_cell_size = handles.max_cell_size;
Diffstack=handles.Diffstack;
img_stacks = handles.img_stacks;
yeast_seg = handles.yeast_seg;
if Diffstack
    images1 = handles.images1;   
    DAPI_images = {images1{handles.indx1}};
    TRANS_images = {images1{handles.indx2}};
    images1 = DAPI_images;
else
    images1 = {handles.sample_file};
    images1
    DAPI_images = {};
    TRANS_images = {};
end
im_prefixes = {'Sample_Seg'};
node_num = handles.node_num;
% [images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                                     %Finds the file names within the fish directory
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
dxy1 = round(sqrt(max_nucleus_size/pi)); 
images1
A1_segment_predefined_variables_streamlined_generalized_Fourth(Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_stacks,yeast_seg,images1,im_prefixes,dxy1,Diffstack,DAPI_images,TRANS_images)
%%% Display segmented boundaries
handles.Xlim_fig = get(handles.Sample_Image_axes,'Xlim');
handles.Ylim_fig = get(handles.Sample_Image_axes,'Ylim');
load([outfile_prefix_seg 'Lab_Sample_Seg'],'cells')
load([outfile_prefix_seg 'nuclei_Sample_Seg'],'Label_mid')
handles.nuc_border = zeros(size(handles.DAPI_ims));
handles.cell_border = zeros(size(handles.DAPI_ims(:,:,1)));
DAPI_ims = handles.DAPI_ims;
TRANS_ims = handles.TRANS_ims;
handles.nuc_border(round(size(DAPI_ims,1)/4):size(DAPI_ims,1)-round(size(DAPI_ims,1)/4),round(size(DAPI_ims,2)/4):size(DAPI_ims,2)-round(size(DAPI_ims,2)/4),:) = Label_mid;
handles.cell_border(round(size(DAPI_ims,1)/4):size(DAPI_ims,1)-round(size(DAPI_ims,1)/4),round(size(DAPI_ims,2)/4):size(DAPI_ims,2)-round(size(DAPI_ims,2)/4)) = cells;
NucBorder1 = bwmorph(handles.nuc_border(:,:,handles.DAPI_slice),'remove')*100000;
NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num);
se = strel('disk',6);
NucBorder1 = imclose(NucBorder1,se); %NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = imclose(NucBorder1,se); NucBorder1 = bwmorph(NucBorder1,'thicken',handles.thicken_num); NucBorder1 = imclose(NucBorder1,se);
CellBorder1 = bwmorph(handles.cell_border,'remove')*100000;
CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num);
se = strel('disk',6);
 CellBorder1 = imclose(CellBorder1,se); 
 if false %handles.thicken_num >1
    CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = imclose(CellBorder1,se); CellBorder1 = bwmorph(CellBorder1,'thicken',handles.thicken_num); CellBorder1 = imclose(CellBorder1,se);
 end
 Mid_DAPI = max(DAPI_ims(:,:,handles.DAPI_slice),[],3);
Nuc1 = Mid_DAPI/max(DAPI_ims(:));
Nuc2 = imadjust(Nuc1,[handles.DAPI_min handles.DAPI_max]);
Trans1 = max(TRANS_ims(:,:,handles.TRANS_slice),[],3);
Trans2 = Trans1/max(TRANS_ims(:));
Trans3 = imadjust(Trans2,[handles.TRANS_min handles.TRANS_max]);
R = Trans3+NucBorder1+CellBorder1; B = Nuc2+Trans3+CellBorder1; G = Trans3+NucBorder1;
RGB = cat(3,R,G,B);
axes(handles.Sample_Image_axes);
imshow(RGB,[]); set(handles.Sample_Image_axes,'Visible','On')
ax = gca;
ax.Toolbar.Visible = 'on';
set(handles.Sample_Image_axes,'Xlim',handles.Xlim_fig)
set(handles.Sample_Image_axes,'Ylim',handles.Ylim_fig)
set(handles.Sample_Image_axes,'xticklabel',[])
set(handles.Sample_Image_axes,'xtick',[])
set(handles.Sample_Image_axes,'yticklabel',[])
set(handles.Sample_Image_axes,'ytick',[])
set(handles.Bound_thick_slider,'SliderStep', [1/11, 0.1])
set(handles.Bound_thick_slider,'Min', 0)
set(handles.Bound_thick_slider,'Max', 10)
set(handles.Bound_thick_slider,'Enable', 'On')
set(handles.Bound_thick_text,'Enable', 'On')
set(handles.Nuc_border_checkbox,'Enable', 'On')
set(handles.Cell_border_checkbox,'Enable', 'On')
set(handles.Nuc_border_checkbox,'Value', 1)
set(handles.Cell_border_checkbox,'Value', 1)
guidata(hObject, handles);


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
thA = handles.thA;
ths = handles.ths;
Yth = handles.Yth;
Ych = handles.Ych;
ManTh = handles.ManTh;
% if ManTh
%     Yth = 1;
% end
Ywin = handles.Ywin;
Yim = handles.Yim;
file_type = handles.file_type;
segment_mode = handles.segment_mode;
Seg_thres = handles.Seg_thres;
fish_dir = handles.fish_dir;
outfile_prefix_seg = handles.outfile_prefix_seg;
node_num = handles.node_num;
max_int_thres = handles.max_int_thres;
max_int_spot_th = handles.max_int_spot_th;
min_nucleus_size = handles.min_nucleus_size;
max_nucleus_size = handles.max_nucleus_size;
min_cell_size = handles.min_cell_size;
max_cell_size = handles.max_cell_size;
img_stacks = handles.img_stacks;
yeast_seg = handles.yeast_seg;
images1 = handles.images1;
im_prefixes = handles.im_prefixes;
node_num = handles.node_num;
Diffstack=handles.Diffstack;
fish_dir = handles.fish_dir;
Diffstack
[images1,folders1] = Get_Dir_All(fish_dir,Ywin);
 clear images2 im_prefixes  
 images2 = {};          %This will have the file names (not the entire directory and not the type of file)
 for i = 1:size(images1,2)
     if Ywin                     %The type of slash depends on whether it is windows
     beginning1 = strfind(images1{i},'\');      
     else
         beginning1 = strfind(images1{i},'/');
     end
     beginning1 = beginning1(end)+1;    %Find the part that is not directory
     if strfind(images1{i},'_MMStack')
         end1 = strfind(images1{i},'_MMStack'); %look for MMStack since it is often at the end of the file names
     elseif strfind(images1{i},'_MMImages')
         end1 = strfind(images1{i},'_MMImages'); %look for MMImages since it is often at the end of the file names
     elseif strfind(images1{i},'.ome')
         end1 = strfind(images1{i},'.ome');        %Look for a period otherwise as the end of the file name
     else
         end1 = strfind(images1{i},'.');        %Look for a period otherwise as the end of the file name
     end
     end1 = end1(1)-1;
     images2{i} = images1{i}(beginning1:end1);  %Populate the cell array with the file names
 end 
 images2';  %For visualizing all files extracted (remove semicolon)
  handles.images1 = images1;
  handles.images2 = images2;
if Diffstack
    clear handles.images_DAPI handles.images_TRANS
    handles.images_DAPI = {};
    handles.images_DAPI = {};
  counter_DAPI = 1;             %counter for marker images found +1 (used to populate both DAPI and TRANS cell arrays)
  counter_TRANS = 0;            %counter for TRANS images found (used to determine if TRANS prefix was correct)
  handles.DAPI_prefix
  for i = 1:size(images2,2) 
      if strfind(images2{i},handles.DAPI_prefix) == 1
          handles.images_DAPI{counter_DAPI} = images1{i};
          im_prefixes{counter_DAPI} =  images2{i}(size(handles.DAPI_prefix,2)+1:end);   %populate imprefixes (used for later segmentation file names)
          for j = 1:size(images2,2)
              if strfind(images2{j},handles.TRANS_prefix) == 1 &...         %Look for TRANS prefix at beginning
                      strcmp(images2{j}(size(handles.TRANS_prefix,2)+1:end),... %Look for the rest of image name to match between TRANS
                      images2{i}(size(handles.DAPI_prefix,2)+1:end))        %and DAPI
                      handles.images_TRANS{counter_DAPI} = images1{j}; 
              end
          end
          counter_TRANS = counter_TRANS+1;                                  %Add to counter (even if match not found)
          counter_DAPI = counter_DAPI+1;                                    %Add to counter (even if match not found)
       if strfind(images2{i},handles.TRANS_prefix) %searching for TRANS prefix
%            handles.images_TRANS{counter_TRANS} = images1{i};
           counter_TRANS = counter_TRANS+1;
       end
      end
  end
  if counter_TRANS == 0
      msgbox('There were no images found with the boundary prefix')
  end
  if isempty(handles.images_DAPI)
      msgbox('There were no images found with the marker prefix')     
  elseif isempty(handles.images_TRANS)
      msgbox('There were no images paired with the prefixes specified')
  end
else
%[images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                %Finds the file names within the fish directory (old way)
im_prefixes = images2;
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
end
dxy1 = round(sqrt(max_nucleus_size/pi));  
DAPI_images = handles.images_DAPI
TRANS_images= handles.images_TRANS
A1_segment_predefined_variables_streamlined_generalized_Fourth(Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_stacks,yeast_seg,images1,im_prefixes,dxy1,Diffstack,DAPI_images,TRANS_images)
msgbox('Segmentation Finished')
% if ManTh
%     Yth = 0;
%     A1_segment_predefined_variables_streamlined_generalized(Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,min_nucleus_size,max_nucleus_size,min_cell_size,max_cell_size,img_stacks,yeast_seg,images1,im_prefixes,dxy1,Diffstack,DAPI_images,TRANS_images)
% end
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Start_RNA_Thresholding_pushbutton.
function Start_RNA_Thresholding_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Start_RNA_Thresholding_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;  %Range of thresholds (e.g. [1,2,100]
ths = handles.ths;
outfile_prefix_RNA = handles.outfile_prefix_seg;
img_stacks = handles.img_stacks;
Yth = handles.Yth;
Ych = handles.Ych;
ManTh = handles.ManTh;
Ywin = handles.Ywin;
Yim = handles.Yim;
file_type = handles.file_type;
segment_mode = handles.segment_mode;
Seg_thres = handles.Seg_thres;
fish_dir = handles.fish_dir;
outfile_prefix_seg = handles.outfile_prefix_seg;
node_num = handles.node_num;
max_int_thres = 1;%handles.max_int_thres; %Set to 1 if you want to use the maximum intensity projection for RNA thresholding visualization
max_int_spot_th = 1;%handles.max_int_spot_th;   % % set to 1 if you want the spot detection to be carried out on the maximum intensity projection during thresholding 
[images1,folders1] = Get_Dir_All(fish_dir,Ywin);
Diffstack = 0;
 clear images2 im_prefixes  
 images2 = {};          %This will have the file names (not the entire directory and not the type of file)
 for i = 1:size(images1,2)
     if Ywin                     %The type of slash depends on whether it is windows
     beginning1 = strfind(images1{i},'\');      
     else
         beginning1 = strfind(images1{i},'/');
     end
     beginning1 = beginning1(end)+1;    %Find the part that is not directory
     if strfind(images1{i},'_MMStack')
         end1 = strfind(images1{i},'_MMStack'); %look for MMStack since it is often at the end of the file names
     elseif strfind(images1{i},'_MMImages')
         end1 = strfind(images1{i},'_MMImages'); %look for MMImages since it is often at the end of the file names
     elseif strfind(images1{i},'.ome')
         end1 = strfind(images1{i},'.ome');        %Look for a period otherwise as the end of the file name
     else
         end1 = strfind(images1{i},'.');        %Look for a period otherwise as the end of the file name
     end
     end1 = end1(end)-1;
     images2{i} = images1{i}(beginning1:end1);  %Populate the cell array with the file names
 end 
 images2';  %For visualizing all files extracted (remove semicolon)
  handles.images1 = images1;
  handles.images2 = images2;
if Diffstack
    clear handles.images_DAPI handles.images_TRANS
    handles.images_DAPI = {};
    handles.images_DAPI = {};
  counter_DAPI = 1;             %counter for marker images found +1 (used to populate both DAPI and TRANS cell arrays)
  counter_TRANS = 0;            %counter for TRANS images found (used to determine if TRANS prefix was correct)
  handles.DAPI_prefix
  for i = 1:size(images2,2) 
      if strfind(images2{i},handles.DAPI_prefix) == 1
          handles.images_DAPI{counter_DAPI} = images1{i};
          im_prefixes{counter_DAPI} =  images2{i}(size(handles.DAPI_prefix,2)+1:end);   %populate imprefixes (used for later segmentation file names)
          for j = 1:size(images2,2)
              if strfind(images2{j},handles.TRANS_prefix) == 1 &...         %Look for TRANS prefix at beginning
                      strcmp(images2{j}(size(handles.TRANS_prefix,2)+1:end),... %Look for the rest of image name to match between TRANS
                      images2{i}(size(handles.DAPI_prefix,2)+1:end))        %and DAPI
                      handles.images_TRANS{counter_DAPI} = images1{j}; 
              end
          end
          counter_TRANS = counter_TRANS+1;                                  %Add to counter (even if match not found)
          counter_DAPI = counter_DAPI+1;                                    %Add to counter (even if match not found)
       if strfind(images2{i},handles.TRANS_prefix) %searching for TRANS prefix
%            handles.images_TRANS{counter_TRANS} = images1{i};
           counter_TRANS = counter_TRANS+1;
       end
      end
  end
  if counter_TRANS == 0
      msgbox('There were no images found with the boundary prefix')
  end
  if isempty(handles.images_DAPI)
      msgbox('There were no images found with the marker prefix')     
  elseif isempty(handles.images_TRANS)
      msgbox('There were no images paired with the prefixes specified')
  end
else
%[images1,im_prefixes] = Get_Dir_ImSt(fish_dir,Ywin)      ;                %Finds the file names within the fish directory (old way)
im_prefixes = images2;
% if Ywin == 0
%     images1 = {images1{node_num+1}};
%     im_prefixes = {im_prefixes{node_num+1}};
% end
end
% Main_RNASpotDetect(img_name, tif_path, out_dir, rna_channel, channel_count, t_min, t_max, true);
Run_RNA_detection_predefined_variables_v2GUI(thA,ths,outfile_prefix_RNA,Yth,Ych,ManTh,Ywin,Yim,file_type,segment_mode,Seg_thres,fish_dir,outfile_prefix_seg,node_num,max_int_thres,max_int_spot_th,handles.images1,im_prefixes,img_stacks)


function Beginning_thresh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Beginning_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;  %Range of thresholds (e.g. [1,2,100]
thA(1) = str2num(get(hObject,'String'))
handles.thA = thA;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Beginning_thresh_edit as text
%        str2double(get(hObject,'String')) returns contents of Beginning_thresh_edit as a double


% --- Executes during object creation, after setting all properties.
function Beginning_thresh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Beginning_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Interval_thresh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Interval_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;  %Range of thresholds (e.g. [1,2,100]
thA(2) = str2num(get(hObject,'String'))
handles.thA = thA;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of Interval_thresh_edit as text
%        str2double(get(hObject,'String')) returns contents of Interval_thresh_edit as a double


% --- Executes during object creation, after setting all properties.
function Interval_thresh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Interval_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function End_thresh_edit_Callback(hObject, eventdata, handles)
% hObject    handle to End_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
thA = handles.thA;  %Range of thresholds (e.g. [1,2,100]
thA(3) = str2num(get(hObject,'String'))
handles.thA = thA;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of End_thresh_edit as text
%        str2double(get(hObject,'String')) returns contents of End_thresh_edit as a double


% --- Executes during object creation, after setting all properties.
function End_thresh_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to End_thresh_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RNA_channel3_edit_Callback(hObject, eventdata, handles)
% hObject    handle to RNA_channel3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ych = handles.Ych;  %Range of thresholds (e.g. [1,2,100]
Ych(6) = str2num(get(hObject,'String'))
handles.Ych = Ych;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of RNA_channel3_edit as text
%        str2double(get(hObject,'String')) returns contents of RNA_channel3_edit as a double


% --- Executes during object creation, after setting all properties.
function RNA_channel3_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RNA_channel3_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RNA_channel2_edit_Callback(hObject, eventdata, handles)
% hObject    handle to RNA_channel2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ych = handles.Ych;  %Range of thresholds (e.g. [1,2,100]
Ych(5) = str2num(get(hObject,'String'))
handles.Ych = Ych;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of RNA_channel2_edit as text
%        str2double(get(hObject,'String')) returns contents of RNA_channel2_edit as a double


% --- Executes during object creation, after setting all properties.
function RNA_channel2_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RNA_channel2_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function RNA_channel1_edit_Callback(hObject, eventdata, handles)
% hObject    handle to RNA_channel1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Ych = handles.Ych;  %Range of thresholds (e.g. [1,2,100]
Ych(4) = str2num(get(hObject,'String'))
handles.Ych = Ych;
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of RNA_channel1_edit as text
%        str2double(get(hObject,'String')) returns contents of RNA_channel1_edit as a double


% --- Executes during object creation, after setting all properties.
function RNA_channel1_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RNA_channel1_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox20.
function checkbox20_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.max_int_thres = get(hObject,'Value');
% Hint: get(hObject,'Value') returns toggle state of checkbox20


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
