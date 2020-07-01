function [tStruct_iw,tStruct_mw,tStruct_mi,clusId,params] = le_mn_clusterPatterns(plotFlag,printFlag)
% This function acts as a wrapper around le_mn_spectral change. 
% 1) Loads ClusScripts across subjects; saves group variable (that can be
% loaded in connectivity script)
% 2) clusters each activation pattern based on features
% 3) Plots on MNI brain

% 12/25 - current set up works nicely. 3 vs 4 k clus. Cluster purely on
% functional activation (t-test w/ instruct, wait and move). It nicely
% separates the activation patterns for HUP084. Need to clean up these
% figures and then plot anatomy

% 12/28 - this is still in development. Has basic funcitonality but needs
% to be significantly cleaned up.

% 1/1/2019 - cleaned up and organized with subfunctions. params struc sets
% config. Need to save figures. moving connectivity analysis to separate
% script (le_mn_connectClusters). see figs saved in printDir from this date
% with associated params file - 5 anatomically and spectrally distributed clusters
% may map onto: salience (MTL/PFC HG), movement (perirolandic broadband),...
% go/no-go (parietofrontal beta), planning?, language?. See how it works
% out across 8 subj. Can consider writing it up +/- connectivity. 

% if both plotFlag == 3 and printFlag == 3, it just returns clusID
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = 1;
end
if ~exist('printFlag','var') || isempty(printFlag) 
    printFlag = 0;
end
%% subject list
% dirs 
dirs = le_dirs;

%% set params
%params.subj_list = {'HUP087'};
params.subj_list = {'HUP084','HUP087','HUP088','HUP089','HUP090','HUP090_1',...
    'HUP091','HUP092','HUP111','HUP112'};
params.saveDir = fullfile(dirs.scratch,'clusStruct');
params.saveFile = ['motornetwork-groupClusStruct-' params.subj_list{:} '.mat'];
params.printDir = fullfile(dirs.scratch,'figs','le_mn_clusterPatterns',date);
params.statsDir = fullfile(dirs.scratch,'stats','le_mn');

% cluster filtering
params.p_thresh_clusStats = 0.001;
params.p_thresh_effSel = 0.01;

% anatomy
params.filt_roi = {'perirolandic','temporal','prefrontal','parietal','MTL','occipital'};
params.filt_hemis = {'','','','','',''};
params.filt_subj = {''};
% params.filt_roi_kol= {'r','b','g','y','o','k'};%

for i =1:length(params.filt_roi)
    params.filt_roiLbl{i}= [params.filt_hemis{i} params.filt_roi{i}];
end


% PCA
params.fig.pca.cLims = [-1 1];

%clusSummary
params.fig.clusSummary.fSize = 16;
params.fig.clusSummary.plotBrainFlag = 0;

% feature selection (sets the features to include in the feature matrix)
params.feat_inclFuncFlag = 1;
params.feat_inclFreqFlag = 0;
params.feat_inclTimeFlag = 0;
params.featLbls = {'zPow move','zPow wait','zPow instruct'};
%params.featLbls = {'zPow move','zPow wait','absDiff(moveRvL)'};
%{'ttest2_instructWait_t','ttest2_moveWait_t','ttest2_moveInstruct_t'};

% k-means
params.k = 5; % if set to empty, it will set automatically based on max jump clus.
params.num_iters = 10000;
params.num_k_to_eval = 10;
params.num_pca_to_use = 5;
% clus anat figure
params.fig.clusAnat.cmap = jet(params.k);
params.fig.clusAnat.markerSize = 100;

% counts
params.counts.foiLbls = {'HFA','gamma','beta','alpha','theta'}; % frequency bands to study for counts analysis

% print
params.fileformat = '-depsc2';
params.resolution = '-r300';%-r0


%% load clusStruct group for all subjects
[tStruct_iw_grp,tStruct_mw_grp,tStruct_mi_grp,...
    clusStruct_iw_grp,clusStruct_mw_grp,clusStruct_mi_grp] = ...
    loadGroupClusStruct_local(params.subj_list,params.saveDir,params.saveFile);
disp('group data LOADED...Ready to Go!!!')

%% retain only electrodes with sig effects
clusStruct_iw_grp = clusStruct_iw_grp([clusStruct_iw_grp.p]<=params.p_thresh_clusStats);
clusStruct_mw_grp = clusStruct_mw_grp([clusStruct_mw_grp.p]<=params.p_thresh_clusStats);
clusStruct_mi_grp = clusStruct_mi_grp([clusStruct_mi_grp.p]<=params.p_thresh_clusStats);


sig_uElbl = [{clusStruct_iw_grp.uElbl},...
    {clusStruct_mw_grp.uElbl},...
    {clusStruct_mi_grp.uElbl}];
sig_uElbl = unique(sig_uElbl);

retIdx = false(1,length({tStruct_iw_grp.uElbl}));
for i = 1:length(sig_uElbl)
    idx = find(strcmp(sig_uElbl{i},{tStruct_iw_grp.uElbl}));
    if ~isempty(idx)
        retIdx(idx) = true;
    end
end

% subselect significant electrodes
tStruct_iw_grp = tStruct_iw_grp(retIdx);
tStruct_mw_grp = tStruct_mw_grp(retIdx);
tStruct_mi_grp = tStruct_mi_grp(retIdx);


%% create specStruct and update params
specStruct.tStruct_iw_grp = tStruct_iw_grp;
specStruct.tStruct_mw_grp = tStruct_mw_grp;
specStruct.tStruct_mi_grp = tStruct_mi_grp;
specStruct.clusStruct_iw_grp = clusStruct_iw_grp;
specStruct.clusStruct_mw_grp = clusStruct_mw_grp;
specStruct.clusStruct_mi_grp = clusStruct_mi_grp;


%% filter specStruct (by region)
[roiStruct] = filterSpecStruct_local(specStruct,params);
clear specStruct;


%% Create Feature Matrix

%global feat mat
featMat = [];
featVec_xtick = []; % this vector can be used to mark transitions in the feature matrix
% make feature matrix of z power in HFA, beta, theta
% rows are electrode
% colums are time in move, instruct, wait
for r = 1:length(params.filt_roiLbl)
    roiLbl = params.filt_roiLbl{r};

    % construct feat Mat
    featStruct.(roiLbl).hfa = [cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zHFA_wait_mean),...
        cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zHFA_instruct_mean),...
        cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zHFA_move_mean)];
    featStruct.(roiLbl).beta = [cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zBeta_wait_mean),...
        cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zBeta_instruct_mean),...
        cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zBeta_move_mean)];
    featStruct.(roiLbl).theta = [cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zTheta_wait_mean),...
        cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zTheta_instruct_mean),...
        cat(1,roiStruct.(roiLbl).tStruct_iw_grp.zTheta_move_mean)];

    featMat = [featMat;[featStruct.(roiLbl).hfa,featStruct.(roiLbl).beta,featStruct.(roiLbl).theta]];
  

    % this vector marks transitions in the feature matrix
    if r == 1
        featVec_xtick = nan(1,size(featMat,2));
        featVec_xtick(size(featStruct.(roiLbl).hfa,2)) = 1;
        featVec_xtick(size(featStruct.(roiLbl).hfa,2)+...
            size(featStruct.(roiLbl).beta,2)) = 1;
    end
    
end


%% do PCA
[coeff,score,latent,~,explained] = pca(featMat);


%% recreate tStruct based on roiStruct to mirror feature matrix indexing
[tStruct_iw,tStruct_mw,tStruct_mi] = roiStruct2tStruct_local(roiStruct,params);

%% k-means
% evaluate num of clusters
jumps = nan(params.num_iters,params.num_k_to_eval);
SUMD = jumps;
parfor i = 1:params.num_k_to_eval
[jumps(i,:),SUMD(i,:)] = JumpClus(score(:,1:params.num_pca_to_use),params.num_k_to_eval);
end
%[jumps,SUMD] = JumpClus(featMat,params.num_k_to_eval);

% choose k if not pre-specified
if isempty(params.k)
    [~,params.k] = max(jumps);
end

% do k-means clustering on this matrix to sort these into
% groups (clusId describes group membership)
[clusId] = kmeans(score(:,1:params.num_pca_to_use),params.k,'replicates',params.num_iters);
%[clusId] = kmeans(featMat,params.k,'replicates',params.num_iters);

%% arrange cluster identification in a reliable manner
arange_param = nan(1,params.k);
for k = 1:params.k
    arange_param(k) = nanmean(score(clusId==k,1),1);
end

[~,sort_k_idx] = sort(arange_param,'descend');
clusId_new = nan(size(clusId));
for k = 1:params.k
    clusId_new(clusId==k) = find(sort_k_idx==k);
end
clusId = clusId_new;
clear clusId_new

%% stop here if you just need to return the cluster ID
if plotFlag ==3 && printFlag == 3% returns cluID
    return
end


%% define instruct and move clusters based on K-means/PCA
instID = 1; % this is the cluster number associated with the instruction group of electrodes
moveID = 3; % this is the cluster number associated with the movement group of electrofes


%% !!!!! OVERWRITE with manual identification of HFA related changes
% isHFA_mw = [];
% isHFA_iw = [];
% isMoveClus = [];
% isInstClus = [];
% for i = 1:length(tStruct_mw)
%     isHFA_mw = [isHFA_mw nanmean(tStruct_mw_grp(i).zHFA_move_mean)];
%     isHFA_iw = [isHFA_iw nanmean(tStruct_iw_grp(i).zHFA_instruct_mean)];
%     
%     %isMoveClus =[isMoveClus;(isHFA_mw)];
%     %isInstClus = [isInstClus;(isHFA_iw&~isHFA_mw)];
% end

% instID = 1; % this is the cluster number associated with the instruction group of electrodes
% moveID = 2; % this is the cluster number associated with the movement group of electrofes
% clusId = zeros(size(clusId));
% clusId(logical(isMoveClus)) = moveID;
% clusId(logical(isInstClus)) = instID;

%% do counts analysis
% do counts analysis a
doCounts_local(tStruct_iw,tStruct_mw,tStruct_mi,clusId,params,instID,moveID);

%% do power analysis
% study power differences between instruction and move clusters  
doPowComp_local(tStruct_iw,tStruct_mw,tStruct_mi,clusId,params,instID,moveID);


%% plot 
if plotFlag
%% plot k-means summary
h = []; % hold all fig handles
[h] = plotK_local(jumps,SUMD,h,params);


%% plot feature mat and PCA sorted
% colorplot of z power
[~,sortIdx] = sort(nanmean(featMat,2));
h(end+1) = swagFig([ 0.1789    0.0575    0.2836    0.8225]);hold all
set(gcf,'name','featureMatrix_unclustered')
imagescnan(featMat(sortIdx,:))
set(gca,'ydir','normal','clim',params.fig.pca.cLims)
colormap jet
swagAxes(gca,16,'Time (ms)','Electrodes','Feature Matrix')
set(gca,'xtick',[],'ytick',[1 size(featMat,1)]);hold all
plot(gca,[find((featVec_xtick==1))' find((featVec_xtick==1))']',...
    repmat(get(gca,'ylim'),[2 1])','--k','linewidth',3)
colorbar

%% feature mat organized by clusters
h(end+1) = swagFig([ 0.1789    0.0575    0.2836    0.8225]);
set(gcf,'name','featureMatrix')
for k = 1:params.k
    subplot(params.k,1,k); hold all
    imagescnan(featMat(clusId==k,:))
    set(gca,'ydir','normal','clim',params.fig.pca.cLims)
    colormap jet
    swagAxes(gca,16,'Time (ms)','Electrodes',['Cluster ' num2str(k)])
    set(gca,'xtick',[],'ytick',[1 sum(clusId==k)])
    plot(gca,[find((featVec_xtick==1))' find((featVec_xtick==1))']',...
    repmat(get(gca,'ylim'),[2 1])','--k','linewidth',3)
    colorbar
end



%%
% scatter plot in PCA space
h(end+1)=swagFig;hold all
set(gcf,'name','PCA')
for i = 1:params.k
    s(i).handle = scatter3(score(clusId == i,1),score(clusId == i,2),score(clusId == i,3),50,'filled');
    s(i).lbl = ['Cluster ' num2str(i)];
end
swagAxes(gca,16,'PCA dimension 1', 'PCA dimension 2');
zlabel('PCA dimension 3');
legend([s.handle],{s.lbl}); 
grid('on')
view(gca,[-15 65])

%% plot cluster summary
[hh] = plotClusterSummary_local(roiStruct,clusId,params,h);
h = [h hh];
%% plot all clusters on brain
[hh] = plotAnat_local(tStruct_iw_grp,clusId,params);
h = [h hh];


%% print group figs
printFigs_local(h,params,plotFlag,printFlag)
    
%% plot examples
x = input('Press y to print example figures, otherwise close all','s');
if strcmp(x,'y')
    fclose('all')
    for i = 1:length(tStruct_iw)
       le_mn_spectralChange(tStruct_iw(i).subj,tStruct_iw(i).eLbl,1,1,[date '-k= ' num2str(clusId(i))]);
    end
else
    fclose('all')
end

end



%% subfunctions

function [] = doPowComp_local(tStruct_iw,tStruct_mw,tStruct_mi,clusId,params,instID,moveID);

% Identify time intervals of interest
tBins = tStruct_iw.tBins;

% retain 
retIdx = (clusId == instID) | (clusId == moveID);
tStruct_iw = tStruct_iw(retIdx);
tStruct_mw = tStruct_mw(retIdx);
tStruct_mi = tStruct_mi(retIdx);
clusId = clusId(retIdx);

% task-related
toi_range = [0 1000];
[~,tIdx_start] = min(abs(tBins-toi_range(1)));
[~,tIdx_end] = min(abs(tBins-toi_range(2)));
tInd = tIdx_start:tIdx_end;

%baseline 
toi_range_bl = [-500 0];
[~,tIdx_start] = min(abs(tBins-toi_range_bl(1)));
[~,tIdx_end] = min(abs(tBins-toi_range_bl(2)));
tInd_bl = tIdx_start:tIdx_end;

% get power matrices for pre-task baseline
HFA_bl = cat(1,tStruct_iw.zHFA_wait_mean); 
HFA_bl = nanmean(HFA_bl(:,tInd_bl),2);

beta_bl = cat(1,tStruct_iw.zBeta_wait_mean);
beta_bl = nanmean(beta_bl(:,tInd_bl),2);

theta_bl = cat(1,tStruct_iw.zTheta_wait_mean);
theta_bl = nanmean(theta_bl(:,tInd_bl),2);

% get power matrices for wait-related power changes
HFA_wait = cat(1,tStruct_iw.zHFA_wait_mean);
HFA_wait = nanmean(HFA_wait(:,tInd),2);

beta_wait = cat(1,tStruct_iw.zBeta_wait_mean);
beta_wait = nanmean(beta_wait(:,tInd),2);

theta_wait = cat(1,tStruct_iw.zTheta_wait_mean);
theta_wait = nanmean(theta_wait(:,tInd),2);


% Get power matrices for instruction-related power changes
HFA_instruct = cat(1,tStruct_iw.zHFA_instruct_mean);
HFA_instruct = nanmean(HFA_instruct(:,tInd),2);

beta_instruct = cat(1,tStruct_iw.zBeta_instruct_mean);
beta_instruct = nanmean(beta_instruct(:,tInd),2);

theta_instruct = cat(1,tStruct_iw.zTheta_instruct_mean);
theta_instruct = nanmean(theta_instruct(:,tInd),2);


% Get power matrices for movement-related power changes
HFA_move = cat(1,tStruct_iw.zHFA_move_mean);
HFA_move = nanmean(HFA_move(:,tInd),2);

beta_move = cat(1,tStruct_iw.zBeta_move_mean);
beta_move = nanmean(beta_move(:,tInd),2);

theta_move = cat(1,tStruct_iw.zTheta_move_mean);
theta_move = nanmean(theta_move(:,tInd),2);

% load pStruct
pStruct = loadPStruct_local(params);
pVec = [];

% Two-factor ANOVA for each freq band
% collect power data for ANOVA
anovDat_hfa = [HFA_bl;HFA_wait;HFA_instruct;HFA_move];
anovDat_beta = [beta_bl;beta_wait;beta_instruct;beta_move];
anovDat_theta = [theta_bl;theta_wait;theta_instruct;theta_move];

%Factor 1: Cluster 1 vs. Cluster 3
fac_clusId = [clusId;clusId;clusId;clusId];

%Factor 2: % baseline vs. wait vs Instruct vs. move,
thisFac = ones(sum(retIdx),1);
fac_task = [thisFac;2*thisFac;3*thisFac;4*thisFac];


% Three-way ANOVA (freq as an additional factor)
thisFac = ones(size(anovDat_hfa));
fac_freq = [thisFac;2*thisFac;3*thisFac];
p_threeWay = anovan([anovDat_hfa;anovDat_beta;anovDat_theta],...
    {fac_freq,[fac_clusId;fac_clusId;fac_clusId],[fac_task;fac_task;fac_task]},...
    'model','full','varnames',{'freq','elecGrp','task'});
pVec = [pVec p_threeWay'];
% Two-Way ANOVA (can do this if needed)
% p_hfa = anovan(anovDat_hfa,{fac_clusId,fac_task},'model','interaction','varnames',{'elecGrp','task'});
% p_beta = anovan(anovDat_beta,{fac_clusId,fac_task},'model','interaction','varnames',{'elecGrp','task'});
% p_theta = anovan(anovDat_theta,{fac_clusId,fac_task},'model','interaction','varnames',{'elecGrp','task'});


% Post-hoc tests
% HFA
% Is instructID instruct HFA > moveID instruct HFA?
[~,p,ci,stats] = ttest2(HFA_instruct(clusId==instID),...
    HFA_instruct(clusId==moveID));
pVec = [pVec p];
% Is instructID move HFA > moveID move HFA?
[~,p,ci,stats] = ttest2(HFA_move(clusId==instID),...
    HFA_move(clusId==moveID));
pVec = [pVec p];

% Beta and theta; is movement-related activity different for instID and
% moveID
[~,p,ci,stats] = ttest2(beta_move(clusId==instID),...
    beta_move(clusId==moveID));
pVec = [pVec p];
[~,p,ci,stats] = ttest2(theta_move(clusId==instID),...
    theta_move(clusId==moveID));
pVec = [pVec p];


% Beta and theta; is instruction-related activity different for instID and
% moveID
[~,p,ci,stats] = ttest2(beta_instruct(clusId==instID),...
    beta_instruct(clusId==moveID));
pVec = [pVec p];
[~,p,ci,stats] = ttest2(theta_instruct(clusId==instID),...
    theta_instruct(clusId==moveID));
pVec = [pVec p];

% Beta and theta; is wait-related activity different for instID and
% moveID
[~,p,ci,stats] = ttest2(beta_wait(clusId==instID),...
    beta_wait(clusId==moveID));
pVec = [pVec p];
[~,p,ci,stats] = ttest2(theta_wait(clusId==instID),...
    theta_wait(clusId==moveID));
pVec = [pVec p];


%  beta; is it increased after the wait cue for instruction and motor
%  sites
% instruction sites
[~,p,ci,stats] = ttest(beta_wait(clusId==instID),...
    beta_bl(clusId==instID));
pVec = [pVec p];
% motor sites
[~,p,ci,stats] = ttest(beta_wait(clusId==moveID),...
    beta_bl(clusId==moveID));
pVec = [pVec p];


% theta: is it increased during wait as compared to the pretask baseline
[~,p,ci,stats] = ttest(theta_wait(clusId==instID),...
    theta_bl(clusId==instID));
pVec = [pVec p];
[~,p,ci,stats] = ttest(theta_wait(clusId==moveID),...
    theta_bl(clusId==moveID));
pVec = [pVec p];


%  beta; is it decreased after the instruct cue for instruction and motor
%  sites
% instruction sites
[~,p,ci,stats] = ttest(beta_instruct(clusId==instID),...
    beta_bl(clusId==instID));
pVec = [pVec p];
% motor sites
[~,p,ci,stats] = ttest(beta_instruct(clusId==moveID),...
    beta_bl(clusId==moveID));
pVec = [pVec p];


% theta: is it increased during instruct as compared to the pretask baseline
[~,p,ci,stats] = ttest(theta_instruct(clusId==instID),...
    theta_bl(clusId==instID));
pVec = [pVec p];
[~,p,ci,stats] = ttest(theta_instruct(clusId==moveID),...
    theta_bl(clusId==moveID));
pVec = [pVec p];



%  beta; is it decreased after the movement cue for instruction and motor
%  sites
% instruction sites
[~,p,ci,stats] = ttest(beta_move(clusId==instID),...
    beta_bl(clusId==instID));
pVec = [pVec p];
% motor sites
[~,p,ci,stats] = ttest(beta_move(clusId==moveID),...
    beta_bl(clusId==moveID));
pVec = [pVec p];


% theta: is it decreased during movement as compared to the pretask baseline
[~,p,ci,stats] = ttest(theta_move(clusId==instID),...
    theta_bl(clusId==instID));
pVec = [pVec p];
[~,p,ci,stats] = ttest(theta_move(clusId==moveID),...
    theta_bl(clusId==moveID));
pVec = [pVec p];


% update pVec
pStruct.anova_clus1vs3 = pVec; 
save('pStruct','pStruct');

function [] = doCounts_local(tStruct_iw,tStruct_mw,tStruct_mi,clusId,params,instID, moveID)
% identify the frequency band and aanatomical labels we want to count
foiLbls = params.counts.foiLbls;
roiLbls = params.filt_roi;
% Initialize variables to hold counts data
countStruct.incMatIW = nan(params.k,length(foiLbls)); % increases during instruction
countStruct.decMatIW = countStruct.incMatIW; % decreases during instruction
countStruct.incMatMW = countStruct.incMatIW; % increases during movement
countStruct.decMatMW = countStruct.incMatIW; % decreases during movement
countStruct.effSelIncIW = nan(params.k,1); % movement-selective increases during instruction
countStruct.effSelDecIW = countStruct.effSelIncIW; % movement-selective decreases during instruction
countStruct.effSelIncMW = countStruct.effSelIncIW; % movement selective increases during movement
countStruct.effSelDecMW = countStruct.effSelIncIW; % movement-selective decreases during movement
countStruct.roi = nan(params.k,length(roiLbls));


% above counters as percentages of number of electrodes in the cluster
countStruct.incMatIW_prct = nan(params.k,length(foiLbls)); % increases during instruction
countStruct.decMatIW_prct = countStruct.incMatIW; % decreases during instruction
countStruct.incMatMW_prct = countStruct.incMatIW; % increases during movement
countStruct.decMatMW_prct = countStruct.incMatIW; % decreases during movement
countStruct.effSelIncIW_prct = nan(params.k,1); % movement-selective increases during instruction
countStruct.effSelDecIW_prct = countStruct.effSelIncIW; % movement-selective decreases during instruction
countStruct.effSelIncMW_prct = countStruct.effSelIncIW; % movement selective increases during movement
countStruct.effSelDecMW_prct = countStruct.effSelIncIW; % movement-selective decreases during movement
countStruct.roi_prct = nan(params.k,length(roiLbls));

% collect all freqBandLbls for largest effect at each electrode
%instruction-related increases
fBandLblsInc_IW = cellfun(@(x) x(1),{tStruct_iw.posClusFreqBandLbl},'UniformOutput',1);
fBandLblsInc_IW(cellfun(@isempty,fBandLblsInc_IW)) = {'NoLbl'};

%instruction-related decreases
fBandLblsDec_IW = cellfun(@(x) x(1),{tStruct_iw.negClusFreqBandLbl},'UniformOutput',1);
fBandLblsDec_IW(cellfun(@isempty,fBandLblsDec_IW)) = {'NoLbl'};

%movement-related increases
fBandLblsInc_MW = cellfun(@(x) x(1),{tStruct_mw.posClusFreqBandLbl},'UniformOutput',1);
fBandLblsInc_MW(cellfun(@isempty,fBandLblsInc_MW)) = {'NoLbl'};

%movement-related decreases
fBandLblsDec_MW = cellfun(@(x) x(1),{tStruct_mw.negClusFreqBandLbl},'UniformOutput',1);
fBandLblsDec_MW(cellfun(@isempty,fBandLblsDec_MW)) = {'NoLbl'};

% collect effector selectivity data for instruction- and movement-related
% increases and decreases
effSelIdxInc_IW = cellfun(@(x) x(1)<=params.p_thresh_effSel,{tStruct_iw.posClusDirSelective_anova_p});
effSelIdxDec_IW = cellfun(@(x) x(1)<=params.p_thresh_effSel,{tStruct_iw.negClusDirSelective_anova_p});
effSelIdxInc_MW = cellfun(@(x) x(1)<=params.p_thresh_effSel,{tStruct_mw.posClusDirSelective_anova_p});
effSelIdxDec_MW = cellfun(@(x) x(1)<=params.p_thresh_effSel,{tStruct_mw.negClusDirSelective_anova_p});


% count total number of electrodes in each cluster and each ROI
n_elec_k = nan(1,length(params.k));
n_elec_roi = nan(1,length(roiLbls));
for k = 1:params.k % looping through each cluster
    
    % count number of electrodes in cluster
    n_elec_k(k) = sum(clusId ==k);
    % counts effector selectivty in each of the four conditions
    countStruct.effSelIncIW(k) = sum(effSelIdxInc_IW(clusId==k)); % increases during instruction
    countStruct.effSelDecIW(k) = sum(effSelIdxDec_IW(clusId==k)); % decreases during instruction
    countStruct.effSelIncMW(k) = sum(effSelIdxInc_MW(clusId==k)); % increases during move
    countStruct.effSelDecMW(k) = sum(effSelIdxDec_MW(clusId==k)); % decreases during move
    
    % percentage of electrodes in cluster w/ effector selectivty in each of the four conditions
    countStruct.effSelIncIW_prct(k) = 100*(sum(effSelIdxInc_IW(clusId==k))/n_elec_k(k)); % increases during instruction
    countStruct.effSelDecIW_prct(k) = 100*(sum(effSelIdxDec_IW(clusId==k))/n_elec_k(k)); % decreases during instruction
    countStruct.effSelIncMW_prct(k) = 100*(sum(effSelIdxInc_MW(clusId==k))/n_elec_k(k)); % increases during move
    countStruct.effSelDecMW_prct(k) = 100*(sum(effSelIdxDec_MW(clusId==k))/n_elec_k(k)); % decreases during move
        
    %looping through frequency bands
    for f = 1:length(foiLbls)
    
        % counts the number of increases during instruction
        countStruct.incMatIW(k,f) = sum(strcmp(fBandLblsInc_IW(clusId==k),foiLbls{f}));
        % counts the number of decreases during instruction
        countStruct.decMatIW(k,f) = sum(strcmp(fBandLblsDec_IW(clusId==k),foiLbls{f}));
        % counts the number of increases during move
        countStruct.incMatMW(k,f) = sum(strcmp(fBandLblsInc_MW(clusId==k),foiLbls{f}));
        % counts the number of decrease during instruction
        countStruct.decMatMW(k,f) = sum(strcmp(fBandLblsDec_MW(clusId==k),foiLbls{f}));
  
        % above values as a percentage of number of electrodes in the
        % cluster
        countStruct.incMatIW_prct(k,f) = 100*(sum(strcmp(fBandLblsInc_IW(clusId==k),foiLbls{f}))/n_elec_k(k));
        countStruct.decMatIW_prct(k,f) = 100*(sum(strcmp(fBandLblsDec_IW(clusId==k),foiLbls{f}))/n_elec_k(k));
        countStruct.incMatMW_prct(k,f) = 100*(sum(strcmp(fBandLblsInc_MW(clusId==k),foiLbls{f}))/n_elec_k(k));
        countStruct.decMatMW_prct(k,f) = 100*(sum(strcmp(fBandLblsDec_MW(clusId==k),foiLbls{f}))/n_elec_k(k));
    
    end
    
    for r = 1:length(roiLbls)
        n_elec_roi(r) = sum(strcmp({tStruct_iw.ROI},roiLbls{r}));
        
        countStruct.roi(k,r) = sum(strcmp({tStruct_iw(clusId==k).ROI},roiLbls{r}));
        
        % This percentage is computed as a percentage of all significant
        % electrodes in the region (not total number of electrodes in the
        % cluster
        countStruct.roi_prct(k,r) = 100*(sum(strcmp({tStruct_iw(clusId==k).ROI},roiLbls{r}))/n_elec_roi(r));
    end
    
end

%
% display tables - Columns are frequency bands, Rows are clusters
%increases during instruction
for k = 1:params.k
    rowNames{k} = ['cluster' num2str(k)];
end
T.incIW = array2table(countStruct.incMatIW,'variablenames',foiLbls,'rownames',rowNames);
T.incIW_prct= array2table(countStruct.incMatIW_prct,'variablenames',foiLbls,'rownames',rowNames);

%decreases during instruction
T.decIW = array2table(countStruct.decMatIW,'variablenames',foiLbls,'rownames',rowNames);
T.decIW_prct = array2table(countStruct.decMatIW_prct,'variablenames',foiLbls,'rownames',rowNames);

%increases during move
T.incMW = array2table(countStruct.incMatMW,'variablenames',foiLbls,'rownames',rowNames);
T.incMW_prct = array2table(countStruct.incMatMW_prct,'variablenames',foiLbls,'rownames',rowNames);

% decreases during move
T.decMW = array2table(countStruct.decMatMW,'variablenames',foiLbls,'rownames',rowNames);
T.decMW_prct = array2table(countStruct.decMatMW_prct,'variablenames',foiLbls,'rownames',rowNames);

%Effector selectivity - Columns are instruction vs. movement, and increase vs. decrease
%  Rows are clusters
colNames = {'InstructIncrease','InstructDecrease',...
    'MovementIncrease','MovementDecrease'};
T.effSel = array2table([countStruct.effSelIncIW countStruct.effSelDecIW ...
    countStruct.effSelIncMW countStruct.effSelDecMW],'variablenames',colNames,'rownames',rowNames);
T.effSel_prct = array2table([countStruct.effSelIncIW_prct countStruct.effSelDecIW_prct ...
    countStruct.effSelIncMW_prct countStruct.effSelDecMW_prct],'variablenames',colNames,'rownames',rowNames);


% anatomy
T.roi = array2table(countStruct.roi, 'variablenames',roiLbls,'rownames',rowNames);
T.roi_prct = array2table(countStruct.roi_prct, 'variablenames',roiLbls,'rownames',rowNames);

%keyboard

% STATS
% load pStruct
pStruct = loadPStruct_local(params);

pVec = [];
 
% PROMINENT FREQ BAND CHANGES ACROSS CLUSTERS
% are HFA increases during instruction non-uniformly distributed across clusters?
% compare to null distribution of cluster counts (n_elec_k)
HFA_inc_o = countStruct.incMatIW(:,1)'+countStruct.incMatMW(:,1)';
[chi2,p,df] = chi2test([HFA_inc_o;n_elec_k]);
pVec = [pVec p];

% are HFA decreases non-uniformly distributed across clusters?
% compare to null distribution of cluster counts (n_elec_k)
HFA_dec_o = countStruct.decMatIW(:,1)' + countStruct.decMatMW(:,1)';
[chi2,p,df] = chi2test([HFA_dec_o;n_elec_k]);
pVec = [pVec p];

% are beta decreases non-uniformly distributed across clusters?
% compare to null distribution of cluster counts (n_elec_k)
beta_dec_o = countStruct.decMatIW(:,3)' + countStruct.decMatMW(:,3)';
[chi2,p,df] = chi2test([beta_dec_o;n_elec_k]);
pVec = [pVec p];

% are theta increases non-uniformly distributed across clusters?
% compare to null distribution of cluster counts (n_elec_k)
theta_inc_o = countStruct.incMatIW(:,1)' + countStruct.incMatMW(:,1)';
[chi2,p,df] = chi2test([theta_inc_o;n_elec_k]);
pVec = [pVec p];



%% instID vs. moveID
% movement selectivity
null_counts = [n_elec_k(instID) n_elec_k(moveID)];
instID_counts_IW= [countStruct.effSelIncIW(instID)+countStruct.effSelDecIW(instID) ...
    sum(clusId==instID)-(countStruct.effSelIncIW(instID)+countStruct.effSelDecIW(instID))];
instID_counts_MW= [countStruct.effSelIncMW(instID)+countStruct.effSelDecMW(instID) ...
    sum(clusId==instID)-(countStruct.effSelIncMW(instID)+countStruct.effSelDecMW(instID))];
moveID_counts_IW= [countStruct.effSelIncIW(moveID)+countStruct.effSelDecIW(moveID) ...
    sum(clusId==moveID)-(countStruct.effSelIncIW(moveID)+countStruct.effSelDecIW(moveID))];
moveID_counts_MW= [countStruct.effSelIncMW(moveID)+countStruct.effSelDecMW(moveID) ...
    sum(clusId==moveID)-(countStruct.effSelIncMW(moveID)+countStruct.effSelDecMW(moveID))];

% are instID and moveID electrodes different in instruction-related movement selectivity? 
[chi2,p,df] = chi2test([instID_counts_IW;moveID_counts_IW]);
pVec = [pVec p];

% are instID and moveID electrodes different in movement-related movement selectivity? 
[chi2,p,df] = chi2test([instID_counts_MW;moveID_counts_MW]);
pVec = [pVec p];

% REGIONAL SELECTIVITY
null_counts = n_elec_roi;
instID_roi = countStruct.roi(instID,:);
moveID_roi = countStruct.roi(moveID,:);
% are instID and moveID different in regional distribution from overall counts?
[chi2,p,df] = chi2test([instID_roi; null_counts]); pVec = [pVec p];
[chi2,p,df] = chi2test([moveID_roi; null_counts]); pVec = [pVec p];

% direct comparison
[chi2,p,df] = chi2test([moveID_roi(:,[1:4,6]); instID_roi(:,[1:4,6])]); pVec = [pVec p];


%% Save pStruct
pStruct.clusCounts = pVec;
save('pStruct','pStruct');


function[tStruct_iw_grp,tStruct_mw_grp,tStruct_mi_grp,...
    clusStruct_iw_grp,clusStruct_mw_grp,clusStruct_mi_grp]...
    = loadGroupClusStruct_local(subj_list,saveDir,saveFile)
% check if group variable is saved
cd_mkdir(saveDir);
skipLoop = 0;
if exist(saveFile,'file')
   load(saveFile) 
   if s == length(subj_list)
       skipLoop = 1;
   end
end
if skipLoop == 0 % we're running the loop
    if ~exist('tStruct_iw_grp','var')
        
        tStruct_iw_grp = struct;
        tStruct_mw_grp = struct;
        tStruct_mi_grp =  struct;
        
    end
    
    if ~exist('s','var')
        s = 1;
    end
    
    for s = s:length(subj_list) % start from where we left off in terms of subject number
        anatStruct = le_loadAnatStruct(subj_list{s});
        for i = 1:length(anatStruct);
           
        [tStruct_iw,tStruct_mw,tStruct_mi,~,~,~,...
            clusStruct_iw,clusStruct_mw,clusStruct_mi]...
            =  le_mn_spectralChange(subj_list{s},anatStruct(i).eLbl,0);            
            
            if (i == 1 && s == 1)||~isfield(tStruct_iw_grp,'subj')
                
                tStruct_iw_grp = tStruct_iw;
                tStruct_mw_grp = tStruct_mw;
                tStruct_mi_grp =  tStruct_mi;
                 
               
            elseif isfield(tStruct_iw_grp,'subj')                   
                    tStruct_iw_grp = [tStruct_iw_grp tStruct_iw];
                    tStruct_mw_grp = [tStruct_mw_grp tStruct_mw];
                    tStruct_mi_grp =  [tStruct_mi_grp tStruct_mi];
                
            end
            if (i == 1 && s==1) || ~isfield(clusStruct_iw_grp,'subj')
                clusStruct_iw_grp = clusStruct_iw;
            elseif isfield(clusStruct_iw,'subj')
                clusStruct_iw_grp = [clusStruct_iw_grp clusStruct_iw];
            end 
            if (i == 1 && s==1) || ~isfield(clusStruct_mw_grp,'subj')
                clusStruct_mw_grp = clusStruct_mw;
            elseif isfield(clusStruct_mw,'subj') 
                clusStruct_mw_grp = [clusStruct_mw_grp clusStruct_mw];
            end
             
             if (i == 1 && s==1) || ~isfield(clusStruct_mi_grp,'subj')
                clusStruct_mi_grp =  clusStruct_mi;
             elseif isfield(clusStruct_mi,'subj')
                clusStruct_mi_grp =  [clusStruct_mi_grp clusStruct_mi];
             end
          disp(['SUBJ ' num2str(s) '.....ELEC' num2str(i) '/' num2str(length(anatStruct))])
        end 
        
        cd_mkdir(saveDir)
        
        save(saveFile,'tStruct_iw_grp','tStruct_mw_grp','tStruct_mi_grp',...
            'clusStruct_iw_grp','clusStruct_mw_grp','clusStruct_mi_grp','s');

            
        disp(['SUBJ ' num2str(s) '....DONE, SAVED']);
    end
else 
    return
end
function[specStruct_ROI,nsubj] = filterSpecStruct_local(specStruct,params)


% anatomy filtering
for i =1:length(params.filt_roi)
    % filter by ROI 
    retIdx = strcmp({specStruct.tStruct_iw_grp.ROI},params.filt_roi{i});
    retIdx_clusIW = strcmp({specStruct.clusStruct_iw_grp.ROI},params.filt_roi{i});
    retIdx_clusMW = strcmp({specStruct.clusStruct_mw_grp.ROI},params.filt_roi{i});
    retIdx_clusMI = strcmp({specStruct.clusStruct_mi_grp.ROI},params.filt_roi{i});
    
    %skip this region if no electrodes in this region
    if sum(retIdx) == 0
        continue
    end
    
    % filter by hemisphere if specified
    if ~isempty(params.filt_hemis{i})
        retIdx = retIdx & strcmp({specStruct.tStruct_iw_grp.hemis},params.filt_hemis{i});
        retIdx_clusIW = retIdx_clusIW & strcmp({specStruct.clusStruct_iw_grp.hemis},params.filt_hemis{i});
        retIdx_clusMW = retIdx_clusMW & strcmp({specStruct.clusStruct_mw_grp.hemis},params.filt_hemis{i});
        retIdx_clusMI = retIdx_clusMI & strcmp({specStruct.clusStruct_mi_grp.hemis},params.filt_hemis{i});
    end
    
    specStruct_ROI.(params.filt_roiLbl{i}).tStruct_iw_grp = specStruct.tStruct_iw_grp(retIdx);
    specStruct_ROI.(params.filt_roiLbl{i}).tStruct_mw_grp = specStruct.tStruct_mw_grp(retIdx);
    specStruct_ROI.(params.filt_roiLbl{i}).tStruct_mi_grp = specStruct.tStruct_mi_grp(retIdx);
    
    nsubj.(params.filt_roiLbl{i}) =  length(unique({specStruct.tStruct_iw_grp(retIdx).subj}));

    % replicate for each clusStruct
    specStruct_ROI.(params.filt_roiLbl{i}).clusStruct_iw_grp = ...
        specStruct.clusStruct_iw_grp(retIdx_clusIW);
    specStruct_ROI.(params.filt_roiLbl{i}).clusStruct_mw_grp = ...
        specStruct.clusStruct_mw_grp(retIdx_clusMW);
    specStruct_ROI.(params.filt_roiLbl{i}).clusStruct_mi_grp = ...
        specStruct.clusStruct_mi_grp(retIdx_clusMI);
end
specStruct_ROI.all = specStruct;
nsubj.all = length(unique({specStruct.tStruct_iw_grp.subj}));

function[h] = plotAnat_local(tStruct_featMat,clusId,params)
% tStruct_iw_grp is the tStruct that has been organized in the same electrode order
% as a feature matrix 
cmap = jet(params.k);

% initialize variables
rightIdx = [tStruct_featMat.x]'>0;
ventralIdx = [strcmp({tStruct_featMat.ROI},'MTL')|strcmp({tStruct_featMat.ROI},'temporal')]';
mni_coords = [cat(1,tStruct_featMat.x) cat(1,tStruct_featMat.y)  cat(1,tStruct_featMat.z)];

% snap electrodes to brain
[mni_coords] = snapCoordsOnBrain_local(mni_coords);

for c = 1:params.k
    % Anatomy
    h(c)=swagFig([ 0.0648    0.3188    0.8000    0.4775]);hold all
    set(gcf,'name',['ClusterAnatomy K = ' num2str(c)])
    
    % identify electrodes in this cluster
    retIdx = clusId == c;
    
    % plot brain and electrodes in subplots
    ax.left = subplot(1,3,1);hold all
    plotBrain_local(.5,[-100 6]);
    scatter3(mni_coords((retIdx&~rightIdx),1),...
    mni_coords((retIdx&~rightIdx),2),...
    mni_coords((retIdx&~rightIdx),3),...
    params.fig.clusAnat.markerSize,'filled','markerfacecolor',cmap(c,:));
    camzoom(1.8)
    
    ax.right = subplot(1,3,2);hold all
    plotBrain_local(.5,[100 6]);    
    scatter3(mni_coords((retIdx&rightIdx),1),...
        mni_coords((retIdx&rightIdx),2),...
        mni_coords((retIdx&rightIdx),3),...
        params.fig.clusAnat.markerSize,'filled','markerfacecolor',cmap(c,:));
    camzoom(1.8)

    
    ax.ventral = subplot(1,3,3);hold all
    plotBrain_local(.5,[180 -90]);
    scatter3(mni_coords((retIdx&ventralIdx),1),...
    mni_coords((retIdx&ventralIdx),2),...
    mni_coords((retIdx&ventralIdx),3),...
    params.fig.clusAnat.markerSize,'filled','markerfacecolor',cmap(c,:));
    camzoom(1.5)
end
function [tStruct_iw,tStruct_mw,tStruct_mi] = roiStruct2tStruct_local(roiStruct,params);

for i = 1:length(params.filt_roiLbl)
    roi = params.filt_roiLbl{i};
    if i == 1
        tStruct_iw = [roiStruct.(roi).tStruct_iw_grp];
        tStruct_mw = [roiStruct.(roi).tStruct_mw_grp];
        tStruct_mi = [roiStruct.(roi).tStruct_mi_grp]; 
    else
        tStruct_iw = [tStruct_iw,roiStruct.(roi).tStruct_iw_grp];
        tStruct_mw = [tStruct_mw,roiStruct.(roi).tStruct_mw_grp];
        tStruct_mi = [tStruct_mi,roiStruct.(roi).tStruct_mi_grp];
    end
end
   
function[h] = plotClusterSummary_local(roiStruct,clusId,params,h)

[tStruct_iw,tStruct_mw,tStruct_mi] = roiStruct2tStruct_local(roiStruct,params);


tBins = tStruct_iw.tBins;
fBins = tStruct_iw.fBins;
fSize = params.fig.clusSummary.fSize;
for k = 1:length(unique(clusId))
    % summary of clusters (functional patterns)
    h(k) = swagFig([0    0.0887    0.9570    0.7538]); hold all
%    set(gcf,'name',['ClusterSummaryK=' num2str(k)])

    set(gcf,'name',['ClusterSummary K = ' num2str(k) ...
        ', n elec = ' num2str(sum(clusId==k)) ', n subj = ' ...
        num2str(length(unique({tStruct_iw(clusId==k).subj}))) ')'])
    
    % tstat IW
    subplot(2,4,1); hold all
    meanT = nanmean(cat(3,tStruct_iw(clusId==k).tMat),3);
    thisH1 = imagescnan(tBins,1:length(fBins),meanT); 
    set(gca,'ydir','normal','clim',[-2 2]);colormap('jet')
    setFreqYTicks(gca,fBins)
    swagAxes(gca,fSize,'','',['t(' tStruct_iw(1).retLbl1 '-' tStruct_iw(1).retLbl2 ')'])
    % plot line
    plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])


    % tstat MW
    subplot(2,4,2); hold all
    meanT = nanmean(cat(3,tStruct_mw(clusId==k).tMat),3);
    thisH1 = imagescnan(tBins,1:length(fBins),meanT); 
    set(gca,'ydir','normal','clim',[-2 2]);colormap('jet')
    setFreqYTicks(gca,fBins)
    swagAxes(gca,fSize,'','',['t(' tStruct_mw(1).retLbl1 '-' tStruct_mw(1).retLbl2 ')'])
    % plot line
    plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])

    % tstat MI
    subplot(2,4,3); hold all
    meanT = nanmean(cat(3,tStruct_mi(clusId==k).tMat),3);
    thisH1 = imagescnan(tBins,1:length(fBins),meanT); 
    set(gca,'ydir','normal','clim',[-2 2]);colormap('jet')
    setFreqYTicks(gca,fBins)
    swagAxes(gca,fSize,'','',['t(' tStruct_mi(1).retLbl1 '-' tStruct_mi(1).retLbl2 ')'])
    % plot line
    plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
 
    % z pow wait
    ax(1) = subplot(2,4,5); hold all
    b(1) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zHFA_wait_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zHFA_wait_mean)),[1 .3 .3]);
    b(2) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zBeta_wait_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zBeta_wait_mean)),[.3 1 .3]);
    b(3) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zTheta_wait_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zTheta_wait_mean)),[.3 .3 1]);
    swagAxes(gca,fSize,'Time from cue (ms)','z power','Wait')
    yL = get(gca,'ylim');
    
    
    % z pow instruct
    ax(2) = subplot(2,4,6); hold all
    b(4) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zHFA_instruct_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zHFA_instruct_mean)),[1 .3 .3]);
    b(5) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zBeta_instruct_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zBeta_instruct_mean)),[.3 1 .3]);
    b(6) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zTheta_instruct_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zTheta_instruct_mean)),[.3 .3 1]);
    swagAxes(gca,fSize,'','','Instruction')
    yL = [yL get(gca,'ylim')];
    
    
    % z pow move
    ax(3) = subplot(2,4,7); hold all
    b(7) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zHFA_move_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zHFA_move_mean)),[1 .3 .3]);
    b(8) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zBeta_move_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zBeta_move_mean)),[.3 1 .3]);
    b(9) = plotSEM(tBins,nanmean(cat(1,tStruct_iw(clusId==k).zTheta_move_mean),1),...
            SEM(cat(1,tStruct_iw(clusId==k).zTheta_move_mean)),[.3 .3 1]);
    swagAxes(gca,fSize,'','','Move')
    yL = [yL get(gca,'ylim')];
    
    for bb = 1:length(b)
        set(b(bb),'facealpha',.75)
    end
    
    for j = 1:length(ax)
        axes(ax(j))
        set(gca,'ylim',[min(yL) max(yL)]);
        plot([0 0],get(gca,'ylim'),'linestyle','--','linewidth',3,'color', [.2 .2 .2])

    end
    
    % for the right most subplots - either print brain plots
    % or bar plots showing anatomical distributoin and effector selectivity
    %or summary
    if params.fig.clusSummary.plotBrainFlag == 1
        % anat
        isRightSide = [tStruct_iw.x]'>0;
        % Right side
        subplot(2,4,4); hold all
        plotBrain_local(.5,[100 6]);
        scatter3([tStruct_iw(isRightSide&clusId==k).x]',[tStruct_iw(isRightSide&clusId==k).y]',[tStruct_iw(isRightSide&clusId==k).z]',100,'filled','r')
        camzoom(2.25)
        % left side
        subplot(2,4,8); hold all
        plotBrain_local(.5,[-100 6]);
        scatter3([tStruct_iw(~isRightSide&clusId==k).x]',[tStruct_iw(~isRightSide&clusId==k).y]',[tStruct_iw(~isRightSide&clusId==k).z]',100,'filled','r')
        %scatter3(-100*ones(sum(~isRightSide&clusId==k),1),[tStruct_iw(~isRightSide&clusId==k).y]',[tStruct_iw(~isRightSide&clusId==k).z]',100,'filled','r')
        camzoom(2.25)
    else
        % percentage of electrodes each region that fall into this cluster
        subplot(2,4,4); hold all
        counts_tot = nan(1,length(params.filt_roi));
        counts_k = nan(1,length(params.filt_roi));
        for r = 1:length(params.filt_roi)
            roi = params.filt_roi{r};
            counts_tot(r) = sum(strcmp(roi,{tStruct_iw.ROI}));
            counts_k(r) = sum(strcmp(roi,{tStruct_iw(clusId==k).ROI}));
        end
        counts_k_prct = (counts_k./counts_tot)*100;
        b = bar(counts_k_prct,'linewidth',1.5);
        set(b,'facecolor',[.7 .7 .7],'edgecolor',[.2 .2 .2],'linewidth',2)
        swagAxes(gca,params.fig.clusSummary.fSize,'',...
            '% of sig. electrodes in each region','');
        set(gca,'ylim',[0 100],'xlim',[.5 length(params.filt_roi)+.5],...
            'xtick',1:length(params.filt_roi),'xticklabel',params.filt_roi,'xticklabelrotation',90)
    
                    
        % percentage of electrodes that show effector selectivity
        subplot(2,4,8); hold all
        counts_thisK_tot = length(tStruct_iw(clusId==k));
        %isEffec selective (incorporate this into tStruct)
        
        [counts_thisK_effSelStruct.posClusIW] = sum(findSelective_local({tStruct_iw(clusId==k).posClusDirSelective_anova_p},params));
        [counts_thisK_effSelStruct.negClusIW] = sum(findSelective_local({tStruct_iw(clusId==k).negClusDirSelective_anova_p},params));

        [counts_thisK_effSelStruct.posClusMW] = sum(findSelective_local({tStruct_mw(clusId==k).posClusDirSelective_anova_p},params));
        [counts_thisK_effSelStruct.negClusMW] = sum(findSelective_local({tStruct_mw(clusId==k).negClusDirSelective_anova_p},params));

        %[counts_thisK_effSelStruct.posClusMI] = sum(findSelective_local({tStruct_mi(clusId==k).posClusDirSelective_anova_p},params));
        %[counts_thisK_effSelStruct.negClusMI] = sum(findSelective_local({tStruct_mi(clusId==k).negClusDirSelective_anova_p},params));
          
        % recreate mat
        counts_thisK_effSel = struct2array(counts_thisK_effSelStruct);
        counts_thisK_effSel_lbls = fieldnames(counts_thisK_effSelStruct);
        
        counts_thiK_effSelprct = (counts_thisK_effSel./counts_thisK_tot)*100;

        b = bar(counts_thiK_effSelprct,'linewidth',1.5);
        set(b,'facecolor',[.7 .7 .7],'edgecolor',[.2 .2 .2],'linewidth',2)
        swagAxes(gca,params.fig.clusSummary.fSize,'',...
            '% movement selective ','');
        set(gca,'ylim',[0 60],'xlim',[.5 length(counts_thisK_effSel_lbls)+.5],...
            'xtick',1:length(counts_thisK_effSel_lbls),'xticklabel',counts_thisK_effSel_lbls,'xticklabelrotation',90)

    end
    
    
end
function[isSel] = findSelective_local(thisSelCell,params)

%thisSelCell = {tStruct_iw.posClusDirSelective_anova_p};
pThresh = params.p_thresh_effSel;

% recreate cell array with only first item
thisSelCell = cellfun(@(x) x(1),thisSelCell);

% 
isSel = thisSelCell<=pThresh;
function[h] = plotK_local(jumps,SUMD,h,params);   
% Evaluate k-means
h = [h swagFig([0.1789    0.5700    0.6547    0.3100])];
set(gcf,'name','kMeansSummary')
subplot(1,2,1)
plot(1:params.num_k_to_eval,nanmean(jumps,1),'-k','linewidth',3);
swagAxes(gca,14,'Number of clusters','Reduction in model error (a.u.)')

subplot(1,2,2)
plot(1:params.num_k_to_eval,nanmean(SUMD,1),'-k','linewidth',3);
swagAxes(gca,14,'Number of clusters','Model error (a.u.)')



function[] = printFigs_local(h,params,plotFlag,printFlag)

% print flag
if plotFlag ==1
    if printFlag
        cd_mkdir(params.printDir)
        save('params','params')
        for i = 1:length(h)
            if i == 12 % this is a hack because it was crashing on fig 13
                keyboard
            end
            figure(h(i));
            print(h(i),get(h(i),'Name'),params.fileformat,params.resolution);
            pause(5)
            close
        end
        %close all
    end
end


% Plot Y ticks for TF plots
function [] = setFreqYTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'ytick',yt,...
'yticklabel',round(fBins(yt)))

function [] = plotBrain_local(facealpha,viewAngle)

if ~exist('facealpha','var') || isempty(facealpha)
    facealpha = 0.2;
end

if ~exist('viewAngle','var') || isempty(viewAngle)
    viewAngle = [90 0];
end

%h=swagFig ([0.1789    0.1738    0.6797    0.7063]);hold all
%set(gcf,'name','Cluster Anatomy')

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
res.hLight = camlight('right');
set(res.hLight,'Color',[.5 .5 .5],'Style','infinite');
lighting gouraud %phong
set(gca,'visible','off')

function [mni_coords] = snapCoordsOnBrain_local(mni_coords)

% snap x
mni_coords(abs(mni_coords(:,1))>70,1) = 70;
% snap y
mni_coords(mni_coords(:,2)>60,2) = 56;
% snap z
mni_coords(mni_coords(:,3)>58,3) = 58;

function [pStruct] = loadPStruct_local(params)

cd_mkdir(params.statsDir)
if exist('pStruct.mat','file')
    load('pStruct.mat')  
else
    pStruct = struct;
end
