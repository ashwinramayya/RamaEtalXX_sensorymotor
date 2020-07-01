function [] = le_mn_plotSubjElecs(printFlag)
% this function plots electrode coverage for each subject

%% initalize
if ~exist('printFlag','var') || isempty(printFlag)
    printFlag = 0;
end

dirs = le_dirs;

params.subj_list = {'HUP084','HUP087','HUP088','HUP089','HUP090','HUP090_1',...
    'HUP091','HUP092','HUP111','HUP112'};

params.filt_roi = {'occipital','temporal','MTL','parietal','perirolandic','prefrontal'};
params.filt_roi_kol = jet(length(params.filt_roi));

%params.filt_roi = {'prefrontal','perirolandic','parietal','temporal','MTL'};
params.printDir = fullfile(dirs.scratch,'figs','le_mn_subjElecs',date);
params.fileformat = '-depsc2';
params.resolution = '-r300';%-r0

%% get anat data
for s = 1:length(params.subj_list)
    subj = params.subj_list{s};
    
    if s == 1
        anatStruct_group = le_loadAnatStruct(subj);
    end
    
    % load anatomical structure
    anatStruct = le_loadAnatStruct(subj);
    anatStruct_group = [anatStruct_group, anatStruct];
    
    disp(s)
end

%% initialize vars
h.right =  plotBrain_local(.5,[100 6]);
h.left = plotBrain_local(.5,[-100 6]);
h.ventral = plotBrain_local(.5,[180 -90]);
rightIdx = [anatStruct_group.x]>0;
ventralIdx = strcmp({anatStruct_group.ROI},'MTL')|strcmp({anatStruct_group.ROI},'temporal');
mni_coords = [cat(1,anatStruct_group.x) cat(1,anatStruct_group.y)  cat(1,anatStruct_group.z)];

%% snap electrodes
[mni_coords] = snapCoordsOnBrain_local(mni_coords);

%% plot electrodes
for i = 1:length(params.filt_roi)
    retIdx = strcmp({anatStruct_group.ROI},params.filt_roi{i});
    numsubj_lbl{i} = [params.filt_roi{i} ' (#subj,#elec) = (' num2str(length(unique({anatStruct_group(retIdx).subj}))) ', ' num2str(sum(retIdx)) ')'];
   
    % plot elecs on right
    figure(h.right)
    s_right(i) = scatter3(mni_coords((retIdx&rightIdx),1),mni_coords((retIdx&rightIdx),2),mni_coords((retIdx&rightIdx),3),75,'filled','markerfacecolor',params.filt_roi_kol(i,:));

    % plot elecs on left
    figure(h.left)
    s_left(i) = scatter3(mni_coords((retIdx&~rightIdx),1),mni_coords((retIdx&~rightIdx),2),mni_coords((retIdx&~rightIdx),3),75,'filled','markerfacecolor',params.filt_roi_kol(i,:));
    
    % plot MTL electrodes
    figure(h.ventral)
    s_ventral(i) = scatter3(mni_coords((retIdx&ventralIdx),1),mni_coords((retIdx&ventralIdx),2),mni_coords((retIdx&ventralIdx),3),75,'filled','markerfacecolor',params.filt_roi_kol(i,:));
    
end

%% set legends
fSize = 12;
figure(h.right);legend(s_right,numsubj_lbl,'fontSize',fSize)
figure(h.left);legend(s_left,numsubj_lbl,'fontSize',fSize)
figure(h.ventral);legend(s_ventral,numsubj_lbl,'fontSize',fSize)

%% print figures
if printFlag == 1
    cd_mkdir(params.printDir)
    print(h.right,'subjElecs_right',params.fileformat,params.resolution);
    print(h.left,'subjElecs_left',params.fileformat,params.resolution);
    print(h.ventral,'subjElecs_ventral',params.fileformat,params.resolution);
    close all
end

function [mni_coords] = snapCoordsOnBrain_local(mni_coords)

% snap x
mni_coords(abs(mni_coords(:,1))>70,1) = 70;
% snap y
mni_coords(mni_coords(:,2)>60,2) = 56;
% snap z
mni_coords(mni_coords(:,3)>58,3) = 58;

function [h] = plotBrain_local(facealpha,viewAngle)

if ~exist('facealpha','var') || isempty(facealpha)
    facealpha = 0.2;
end

if ~exist('viewAngle','var') || isempty(viewAngle)
    viewAngle = [90 0];
end

h=swagFig ([0.1789    0.1738    0.6797    0.7063]);hold all
set(gcf,'name','Cluster Anatomy')

% read in the faces and vertices of the surface
% get path to pictures
picpath = fileparts(which('tal3d'));
load(fullfile(fileparts(which('tal3d')),'mni_cortical_surface.mat'));
% draw the cortical surface patches  
hs = patch('faces',f,'vertices',v,'edgecolor','none','facecolor',[.9 .9 .9],'facealpha',facealpha);
% set aspect, view, and light
daspect([1 1 1]);
view(viewAngle)
% save the handle
res.hBrain = hs;
% set lighting
res.hLight = camlight;
set(res.hLight,'Color',[1 1 1],'Style','infinite');
lighting phong
set(gca,'visible','off')