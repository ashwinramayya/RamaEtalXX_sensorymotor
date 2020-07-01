function [tStruct_iw,tStruct_mw,tStruct_mi,...
    powStruct_iw,powStruct_mw,powStruct_mi,...
    clusStruct_iw,clusStruct_mw,clusStruct_mi] = le_mn_spectralChange(subj,eLbl,plotFlag,printFlag,printSubDir)
%* Goal: Obtain functional activity for each electrode (nodal measures; step 1)
%* Returns a structure with clusters of spectral change from each electrode
%%* Also tag each cluster with subj, sess, and electrode (possibly,  anatomical information)
% Generalize to allow for tilt regression and phase reset
%* Generalize to combine across electrodes for (ROI analysis or network analysis)

% Inputs
% subj  =        'HUP084'
% eLbl  = '9-10'; % don't use multiple electrode list as an input -
% 12/23/18 AGR

% load dirs and config
dirs = le_dirs;

% parse inputs
if ~exist('printFolderLabel','var') || isempty(printFolderLabel)
    printFolderLabel = '';
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = 1;
end
if ~exist('printFlag','var') || isempty(printFlag)
    printFlag = 0;
end
if ~exist('printSubDir','var') || isempty(printSubDir)
    printSubDir = '';
end
if ~exist('eLbl','var') || isempty(eLbl)
    eLbl = {'1-2'};
end

% get config
config_mn = getConfig_local('all');

% load anat struct
anatStruct = le_loadAnatStruct(subj,eLbl);

% get power with cluster stats
[tStruct_iw,tStruct_mw,tStruct_mi,config_pow_all,h,thisAx] = getTstruct_local(subj,eLbl,config_mn,plotFlag,anatStruct);

% get general power-related feaures for this electrode
% update tStruct with important data for power analyses
[powStruct_iw,tStruct_iw,h,thisAx] = getPowStruct_local(tStruct_iw,anatStruct,config_pow_all,config_mn,h,thisAx,plotFlag);
[powStruct_mw,tStruct_mw,h,thisAx] = getPowStruct_local(tStruct_mw,anatStruct,config_pow_all,config_mn,h,thisAx,plotFlag);
[powStruct_mi,tStruct_mi,h,thisAx] = getPowStruct_local(tStruct_mi,anatStruct,config_pow_all,config_mn,h,thisAx,plotFlag);


% get features of each sig cluster (move-wait) and update tStruct
% such that we categorize clusters based on freq band and time band
[clusStruct_iw,tStruct_iw] = getSpecFeatures_local(tStruct_iw,powStruct_iw,anatStruct,config_mn,config_pow_all,plotFlag,thisAx.iw,h.iw,thisAx.tf,h.tf);
[clusStruct_mw,tStruct_mw] = getSpecFeatures_local(tStruct_mw,powStruct_mw,anatStruct,config_mn,config_pow_all,plotFlag,thisAx.mw,h.mw,thisAx.tf,h.tf);
[clusStruct_mi,tStruct_mi] = getSpecFeatures_local(tStruct_mi,powStruct_mi,anatStruct,config_mn,config_pow_all,plotFlag,thisAx.mi,h.mi,thisAx.tf,h.tf);

% print flag
if plotFlag ==1
    
    %update TF plot: - plot anatomy, place color bar and labels with anatomy
    updateTFplot_local(thisAx.tf,config_mn.clim_tfDiff,anatStruct)

    if printFlag
        printSubDir = fullfile(dirs.scratch,'figs','le_mn_spectralChange',printSubDir,subj);
        cd_mkdir(printSubDir)
        figLbls = fieldnames(h);
        
        if ~isempty(printSubDir) % if print SubDir, print only TF
            print(h.tf,get(gcf,'Name'),'-depsc2','-r300');
        else
            for i = 1:length(figLbls)
                figure(h.(figLbls{i}));
                print(gcf,get(gcf,'Name'),'-dpng');
            end
        end
        close all
    end
end

fclose('all')


%% SUBFCN
% local get spectral features 
function [clusStruct,tStruct] = getSpecFeatures_local(tStruct,powStruct,anatStruct,config_mn,config_pow_all,plotFlag,thisAx,thisH,thisAx_tf,h_tf);

clusStruct = struct;

% plot false positive rate of cluster statistics
if plotFlag == 1
    figure(thisH)
    thisAx(8) = subplot(2,4,8); hold all
    h1 = histogram(tStruct.shufStats_pos,'facecolor',[1 0 0]); 
    h2 = histogram(tStruct.shufStats_neg,'facecolor',[0 0 1]); 
end

% identify positive stats
clus_pos_idx = find([tStruct.posClusStatsP]<=config_mn.fdr_p_thresh);

% initialize freqBand and timeBand labels (and anovaDir selectivity)
tStruct.posClusFreqBandLbl = cell(size(tStruct.posClusStatsP));
tStruct.negClusFreqBandLbl = cell(size(tStruct.negClusStatsP));
tStruct.posClusTimeBandLbl = cell(size(tStruct.posClusStatsP));
tStruct.negClusTimeBandLbl = cell(size(tStruct.negClusStatsP));
tStruct.posClusDirSelective_anova_p = nan(size(tStruct.posClusStatsP));
tStruct.negClusDirSelective_anova_p = nan(size(tStruct.negClusStatsP));

if ~isempty(clus_pos_idx)
    for c = clus_pos_idx
        
        % mark this clusters statistic on the histogram
        if plotFlag == 1
            figure(thisH)
            thisAx(8) = subplot(2,4,8);hold all
            plot([tStruct.posClusStats(c) tStruct.posClusStats(c)],get(gca,'ylim'),'r','linewidth',3);
        end
        
        if c == 1 % then plot pow spectra
            [clusStruct,tStruct] = updateClusStruct_local(clusStruct,tStruct,powStruct,anatStruct,'pos',c,config_mn,config_pow_all,plotFlag,thisAx,thisH,thisAx_tf,h_tf);
        else % plotFlag is hard coded to off
            [clusStruct,tStruct] = updateClusStruct_local(clusStruct,tStruct,powStruct,anatStruct,'pos',c,config_mn,config_pow_all,0,thisAx,thisH,thisAx_tf,h_tf);
        end
    end
end
clus_neg_idx = find([tStruct.negClusStatsP]<=config_mn.fdr_p_thresh);
if ~isempty(clus_neg_idx)
    for c = clus_neg_idx
        
        % mark this clusters statistic on the histogram
        if plotFlag == 1
            thisAx(8) = subplot(2,4,8); hold all
            plot([tStruct.negClusStats(c) tStruct.negClusStats(c)],get(gca,'ylim'),'b','linewidth',3);
        end
        
        if c == 1 && isempty(clus_pos_idx)% then plot pow spectra
            [clusStruct,tStruct] = updateClusStruct_local(clusStruct,tStruct,powStruct,anatStruct,'neg',c,config_mn,config_pow_all,plotFlag,thisAx,thisH,thisAx_tf,h_tf);
        else % plotFlag is hard coded to off
            [clusStruct,tStruct] = updateClusStruct_local(clusStruct,tStruct,powStruct,anatStruct,'neg',c,config_mn,config_pow_all,0,thisAx,thisH,thisAx_tf,h_tf);
        end
    end
end

    
function  [clusStruct,tStruct] = updateClusStruct_local(clusStruct,tStruct,powStruct,anatStruct,posNegStr,c,config_mn,config_pow_all,plotFlag,thisAx,thisH,thisAx_tf,h_tf)

flbls = fieldnames(tStruct);
albls = fieldnames(anatStruct);
thisClusStruct = struct();
    
% identify clusStat
for f = 1:length(flbls)
    switch flbls{f}
        case {'subj','sessLbl','eLbl','comparison','task','retLbl1','retLbl2','tMat','tInd','fInd'}
            thisClusStruct.(flbls{f}) = tStruct.(flbls{f}); 
    end
end

% update this clusStruct with anatomical info

for f = 1:length(albls)
    switch albls{f}
        case {'anat','x','y','z','hemis','elecfullname','ROI'}
            thisClusStruct.(albls{f}) = anatStruct.(albls{f});
            %tStruct.(albls{f}) = anatStruct.(albls{f});
    end
end

%  unique identifier
thisClusStruct.uElbl = [thisClusStruct.subj '-' thisClusStruct.eLbl];

% thisClus data
% pos v neg
%tBins = mean(config_pow_all.timeBins,2);tBins = tBins(tStruct(1).tInd);
%fBins = config_pow_all.freqBins;fBins = fBins(tStruct(1).fInd); % filter by freq range used for this tStruct
thisClusStruct.tBins = tStruct.tBins;
thisClusStruct.fBins = tStruct.fBins;

switch  posNegStr
    case 'pos'
        thisClusStruct.clusPolarity = 1;
        thisClusStat = tStruct.posClusStats(c);
        thisClusP = tStruct.posClusStatsP(c);
        thisClusInd = tStruct.posClusInd{c};
    case 'neg'
        thisClusStruct.clusPolarity = -1;
        thisClusStat = tStruct.negClusStats(c);
        thisClusP = tStruct.negClusStatsP(c);
        thisClusInd = tStruct.negClusInd{c};
end

%clus data
thisClusStruct.stat = thisClusStat;
thisClusStruct.p = thisClusP;
if strcmp(thisClusStruct.comparison,'instructWait')
    thisClusStruct.InstructOrMove = 0;
elseif strcmp(thisClusStruct.comparison,'moveWait')
    thisClusStruct.InstructOrMove = 1;
elseif strcmp(thisClusStruct.comparison,'moveInstruct')
    thisClusStruct.InstructOrMove = 2;
end

% freq and time cluster
[thisClusStruct.freqIdx,thisClusStruct.timeIdx] = ...
    ind2sub(size(tStruct.tMat),thisClusInd);

thisClusStruct.freqs = tStruct.fBins(thisClusStruct.freqIdx)';
thisClusStruct.times = tStruct.tBins(thisClusStruct.timeIdx);
thisClusStruct.freqMax = max(thisClusStruct.freqs);
thisClusStruct.freqMin = min(thisClusStruct.freqs);
thisClusStruct.freqSpan = max(thisClusStruct.freqs)-min(thisClusStruct.freqs);
thisClusStruct.timeMax = max(thisClusStruct.times);
thisClusStruct.timeMin = min(thisClusStruct.times);

% peak
[thisClusStruct.clusPeaktStat, peakIdxWithinClus]...
    = max(tStruct.tMat(thisClusInd));

%get freq and time data
[thisClusStruct.peakFIdx,thisClusStruct.peakTIdx] = ind2sub(size(tStruct.tMat),thisClusInd(peakIdxWithinClus));
thisClusStruct.peakFreq = tStruct.fBins(thisClusStruct.peakFIdx);
thisClusStruct.peakTime = tStruct.tBins(thisClusStruct.peakTIdx);

% % add zPowDiff mat and vec
thisClusStruct.zPowDiffVecMean = nanmean(tStruct.zPowDiffMat(thisClusStruct.freqIdx,:),1);
thisClusStruct.zPowDiffClusMean =   nanmean(nanmean(tStruct.zPowDiffMat(thisClusStruct.freqIdx,thisClusStruct.timeIdx),1));


% % assign a frequncy band and time band label and update tStruct as well
fbandLbls = fieldnames(config_mn.freqBandLims);
fRange = [thisClusStruct.freqMin thisClusStruct.freqMax];
freqBandLims_mat = struct2cell(config_mn.freqBandLims);
freqBandLims_mat = cat(1,freqBandLims_mat{:});

tBandLbls = fieldnames(config_mn.timeBandlims);
tRange = [thisClusStruct.timeMin thisClusStruct.timeMax];
timeBandLims_mat = struct2cell(config_mn.timeBandlims);
timeBandLims_mat = cat(1,timeBandLims_mat{:});


%template match to a freq band by calculating distance between
% between cluster's frequency range (span and mean frequency) 
% and each frequency band
% D(i,j) is the distance between the ith observation X (activations),...
% and the jth observation in Y (freq bands)
D = pdist2(fRange,freqBandLims_mat);
[~,fBandIdx] = min(D,[],2);

Dt = pdist2(tRange,timeBandLims_mat);
[~,tBandIdx] = min(Dt,[],2);

%update clus
thisClusStruct.freqBandLbl = fbandLbls{fBandIdx};
thisClusStruct.timeBandLbl = tBandLbls{tBandIdx};

% clusLbls
thisClusStruct.clusLbl = ['isMove' num2str(thisClusStruct.InstructOrMove) '-F' num2str(round(thisClusStruct.peakFreq)) '-T' num2str(round(thisClusStruct.peakTime))];

[thisClusStruct] = getClusSelectivity(thisClusStruct,powStruct,config_pow_all,config_mn,plotFlag,thisAx,thisH,thisAx_tf,h_tf);

if isfield(clusStruct,'subj')        
    clusStruct = [clusStruct thisClusStruct];
else
    clusStruct = thisClusStruct;
end

% update tStruct
switch  posNegStr
    case 'pos'
       tStruct.posClusFreqBandLbl{c} = fbandLbls{fBandIdx};
       tStruct.posClusTimeBandLbl{c} = tBandLbls{tBandIdx};
       tStruct.posClusDirSelective_anova_p(c) = thisClusStruct.anova_dirSelectivity_p;  
       
    case 'neg'
       tStruct.negClusFreqBandLbl{c} = fbandLbls{fBandIdx};
       tStruct.negClusTimeBandLbl{c} = tBandLbls{tBandIdx};
       tStruct.negClusDirSelective_anova_p(c) = thisClusStruct.anova_dirSelectivity_p;  
end

function [powStruct,tStruct,h,thisAx] = getPowStruct_local(tStruct,anatStruct,config_pow_all,config_mn,h,thisAx,plotFlag);
subj = tStruct.subj;
switch tStruct.sessLbl
    case 'all'
        sessList = [];
end
elbl = tStruct.eLbl;
task = tStruct.task;


% get values for x-axis (time bins) and yaxis (freq bins, accounting for ...
%custom fRange used to generate this tStruct)
% tBins = nanmean(config_pow_all.timeBins,2)'; 
% tBins = tBins(tStruct(1).tInd);
% fBins = [config_pow_all.freqBins];
% fBins = fBins(tStruct(1).fInd); % filter by freq range used for this tStruct

% load events
[events] = le_loadEvents(subj,sessList,task);

% get task conditions for instruct, move, wait x (right, left and mouth)
% main effect
retStruct.instruct = strcmp({events.type},'INSTRUCT');
retStruct.wait = strcmp({events.type},'WAIT');
retStruct.move = strcmp({events.type},'MOVE');

% get task conditions for instruct, move, wait x (right, left and mouth)
% main effect
retStruct.instruct = strcmp({events.type},'INSTRUCT');
retStruct.wait = strcmp({events.type},'WAIT');
retStruct.move = strcmp({events.type},'MOVE');

% main effect
retStruct.left = strcmp({events.item},'Left Hand');
retStruct.mouth = strcmp({events.item},'Mouth and Tongue');
retStruct.right = strcmp({events.item},'Right Hand');

% interaction
retStruct.instructleft = strcmp({events.type},'INSTRUCT')&strcmp({events.item},'Left Hand');
retStruct.waitleft = strcmp({events.type},'WAIT')&strcmp({events.item},'Left Hand');
retStruct.moveleft =  strcmp({events.type},'MOVE')&strcmp({events.item},'Left Hand');

retStruct.instructmouth = strcmp({events.type},'INSTRUCT')&strcmp({events.item},'Mouth and Tongue');
retStruct.waitmouth = strcmp({events.type},'WAIT')&strcmp({events.item},'Mouth and Tongue');
retStruct.movemouth =  strcmp({events.type},'MOVE')&strcmp({events.item},'Mouth and Tongue');

retStruct.instructright = strcmp({events.type},'INSTRUCT')&strcmp({events.item},'Right Hand');
retStruct.waitright = strcmp({events.type},'WAIT')&strcmp({events.item},'Right Hand');
retStruct.moveright =  strcmp({events.type},'MOVE')&strcmp({events.item},'Right Hand');

% get Pow for all events and collapse within cluster
[subjPow,~,subjRawPow] = le_calcPow(subj,sessList,elbl,task,[],[]);
retLbls = fieldnames(retStruct);
% filter subjPow by the freq and time index used in the t-matrix
% calculation
subjPow = subjPow(tStruct.fInd,tStruct.tInd,:);
subjRawPow = subjRawPow(tStruct.fInd,tStruct.tInd,:);

% remove inf trials here:
% This is a hack: occurs very rarely, one e.g., HUP090_1
subjPow(isinf(subjPow)) = nan;
subjRawPow(isinf(subjRawPow)) = nan;

% make zPowVec for HFA and LFA
[~,fIdx_HFA(1)] = min(abs(tStruct.fBins - config_mn.freqDef.HFA(1)));
[~,fIdx_HFA(2)] = min(abs(tStruct.fBins - config_mn.freqDef.HFA(2)));
zClusPowVec_HFA = shiftdim(nanmean(subjPow(unique(fIdx_HFA(1):fIdx_HFA(2)),:,:),1),1)'; 
zClusPowVec_HFA = (zClusPowVec_HFA-nanmean(nanmean(zClusPowVec_HFA,2)))/nanstd(nanmean(zClusPowVec_HFA,2));

% beta 
[~,fIdx_beta(1)] = min(abs(tStruct.fBins - config_mn.freqBandLims.beta(1)));
[~,fIdx_beta(2)] = min(abs(tStruct.fBins - config_mn.freqBandLims.beta(2)));
zClusPowVec_beta = shiftdim(nanmean(subjPow(unique(fIdx_beta(1):fIdx_beta(2)),:,:),1),1)'; 
zClusPowVec_beta = (zClusPowVec_beta-nanmean(nanmean(zClusPowVec_beta,2)))/nanstd(nanmean(zClusPowVec_beta,2));


% theta
[~,fIdx_theta(1)] = min(abs(tStruct.fBins - config_mn.freqBandLims.theta(1)));
[~,fIdx_theta(2)] = min(abs(tStruct.fBins - config_mn.freqBandLims.theta(2)));
zClusPowVec_theta = shiftdim(nanmean(subjPow(unique(fIdx_theta(1):fIdx_theta(2)),:,:),1),1)'; 
zClusPowVec_theta = (zClusPowVec_theta-nanmean(nanmean(zClusPowVec_theta,2)))/nanstd(nanmean(zClusPowVec_theta,2));


%LFA
[~,fIdx_LFA(1)] = min(abs(tStruct.fBins - config_mn.freqDef.LFA(1)));
[~,fIdx_LFA(2)] = min(abs(tStruct.fBins - config_mn.freqDef.LFA(2)));
zClusPowVec_LFA = shiftdim(nanmean(subjPow(unique(fIdx_LFA(1):fIdx_LFA(2)),:,:),1),1)'; 
zClusPowVec_LFA = (zClusPowVec_LFA-nanmean(nanmean(zClusPowVec_LFA,2)))/nanstd(nanmean(zClusPowVec_LFA,2));

% add zPowDiff mat and vec
zPowDiffMat = nanmean(tStruct.pow1,3) - nanmean(tStruct.pow2,3);
zPowDiffVec_HFA = nanmean(zPowDiffMat(unique(fIdx_HFA(1):fIdx_HFA(2)),:),1);
zPowDiffVec_LFA = nanmean(zPowDiffMat(unique(fIdx_LFA(1):fIdx_LFA(2)),:),1);

% plot this clus powSpectra below so that we get the whole picture for the
% spectral change for when the t-statistics were computed interval
powSpectra = shiftdim(nanmean(subjRawPow(:,:,:),2),2)';
%% (   ) Implement fooof to parameterize power spectra
%Haller M, Donoghue T, Peterson E, Varma P, Sebastian P, Gao R, Noto T, Knight RT, Shestyuk A,
%Voytek B (2018) Parameterizing Neural Power Spectra. bioRxiv, 299859.
%doi: https://doi.org/10.1101/299859
%paramPowStruct = parametrizePowSpectra_local(clusPowSpectra);
%% (   ) collect phase, EEG,at peak frequency and time, calculate get phase reset
%[subjPhase] = le_calcPhase(subj,sessList,elbl,task,[],[]);
%keyboard

%Get EEG  across all trials
%no buffer 
bufferMS = 0;
%bufferMS = config_pow_all.bufferMS
% shortened time range
switch tStruct.comparison
    case 'instructWait'
        tRange = config_mn.tRange_iw;
    case 'moveWait'
        tRange = config_mn.tRange_mw;
    case 'moveInstruct'
        tRange = config_mn.tRange_mi;
end


durationMS = tRange(end)-tRange(1);
offsetMS = tRange(1);
EEG = gete_ms_wrapper(elbl,events, ...
                   durationMS,offsetMS,bufferMS,...
                   config_pow_all.filtfreq, config_pow_all.filttype, config_pow_all.filtorder);


% update powStruct
powStruct.retStruct = retStruct;
powStruct.retLbls = retLbls;
powStruct.zClusPowVec_HFA = zClusPowVec_HFA;
powStruct.zClusPowVec_LFA = zClusPowVec_LFA;
powStruct.subjPow = subjPow;
powStruct.subjRawPow = subjRawPow;
powStruct.powSpectra = powSpectra;
powStruct.EEG = EEG;

% tStruct
tStruct.zPowDiffMat = zPowDiffMat;
tStruct.zPowDiffVec_HFA = zPowDiffVec_HFA;
tStruct.zPowDiffVec_LFA = zPowDiffVec_LFA;

tStruct.zHFA_wait_mean = nanmean(zClusPowVec_HFA(retStruct.wait,:),1);
tStruct.zHFA_instruct_mean = nanmean(zClusPowVec_HFA(retStruct.instruct,:),1);
tStruct.zHFA_move_mean = nanmean(zClusPowVec_HFA(retStruct.move,:),1);
tStruct.zHFA_wait_sem = SEM(zClusPowVec_HFA(retStruct.wait,:));
tStruct.zHFA_instruct_sem = SEM(zClusPowVec_HFA(retStruct.instruct,:));
tStruct.zHFA_move_sem = SEM(zClusPowVec_HFA(retStruct.move,:));

tStruct.zBeta_wait_mean = nanmean(zClusPowVec_beta(retStruct.wait,:),1);
tStruct.zBeta_instruct_mean = nanmean(zClusPowVec_beta(retStruct.instruct,:),1);
tStruct.zBeta_move_mean = nanmean(zClusPowVec_beta(retStruct.move,:),1);
tStruct.zBeta_wait_sem = SEM(zClusPowVec_beta(retStruct.wait,:));
tStruct.zBeta_instruct_sem = SEM(zClusPowVec_beta(retStruct.instruct,:));
tStruct.zBeta_move_sem = SEM(zClusPowVec_beta(retStruct.move,:));

tStruct.zTheta_wait_mean = nanmean(zClusPowVec_theta(retStruct.wait,:),1);
tStruct.zTheta_instruct_mean = nanmean(zClusPowVec_theta(retStruct.instruct,:),1);
tStruct.zTheta_move_mean = nanmean(zClusPowVec_theta(retStruct.move,:),1);
tStruct.zTheta_wait_sem = SEM(zClusPowVec_theta(retStruct.wait,:));
tStruct.zTheta_instruct_sem = SEM(zClusPowVec_theta(retStruct.instruct,:));
tStruct.zTheta_move_sem = SEM(zClusPowVec_theta(retStruct.move,:));


tStruct.zLFA_wait_mean = nanmean(zClusPowVec_LFA(retStruct.wait,:),1);
tStruct.zLFA_instruct_mean = nanmean(zClusPowVec_LFA(retStruct.instruct,:),1);
tStruct.zLFA_move_mean = nanmean(zClusPowVec_LFA(retStruct.move,:),1);
tStruct.zLFA_wait_sem = SEM(zClusPowVec_LFA(retStruct.wait,:));
tStruct.zLFA_instruct_sem = SEM(zClusPowVec_LFA(retStruct.instruct,:));
tStruct.zLFA_move_sem = SEM(zClusPowVec_LFA(retStruct.move,:));


% update powStruct with anatomical struct
albls = fieldnames(anatStruct);
for f = 1:length(albls)
    switch albls{f}
        case {'anat','x','y','z','hemis','elecfullname','ROI'}
            powStruct.(albls{f}) = anatStruct.(albls{f});
     end
end
%% plot pow spectra
if plotFlag
    fSize = 22;
    
    figure(h.tf);

    % plot HFA/LFA
    switch tStruct.comparison
        case 'instructWait' % this is instruct-wait 
        aLbl = 'iw';
        
        % plot WAIT
        axes(thisAx.tf(1)); hold all    

        % HFA
        b(1)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.waitleft,:),1),...
            SEM(zClusPowVec_HFA(retStruct.waitleft,:)),[1 .3 .3]);
 
        b(2)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.waitright,:),1),...
            SEM(zClusPowVec_HFA(retStruct.waitright,:)),[1 .5 .5]);

        b(3)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.waitmouth,:),1),...
            SEM(zClusPowVec_HFA(retStruct.waitmouth,:)),[1 .7 .7]);

        % LFA
        b(4)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.waitleft,:),1),...
            SEM(zClusPowVec_LFA(retStruct.waitleft,:)),[.3 .3 1]);
 
        b(5)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.waitright,:),1),...
            SEM(zClusPowVec_LFA(retStruct.waitright,:)),[.5 .5 1]);

        b(6)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.waitmouth,:),1),...
            SEM(zClusPowVec_LFA(retStruct.waitmouth,:)),[.7 .7 1]);
      
            
        %swagAxes(gca,fSize,'Time from WAIT cue (ms)',['z power ' tStruct.retLbl2])
        swagAxes(gca,fSize,'Time from cue (ms)','z power','Wait')
  
        for bb = 1:length(b)
            set(b(bb),'facealpha',.75)
        end

        % plot INSTRUCT
        axes(thisAx.tf(2)); hold all    

        % HFA
        b(1)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.instructleft,:),1),...
            SEM(zClusPowVec_HFA(retStruct.instructleft,:)),[1 .3 .3]);
 
        b(2)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.instructright,:),1),...
            SEM(zClusPowVec_HFA(retStruct.instructright,:)),[1 .5 .5]);

        b(3)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.instructmouth,:),1),...
            SEM(zClusPowVec_HFA(retStruct.instructmouth,:)),[1 .7 .7]);

        % LFA
        b(4)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.instructleft,:),1),...
            SEM(zClusPowVec_LFA(retStruct.instructleft,:)),[.3 .3 1]);
 
        b(5)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.instructright,:),1),...
            SEM(zClusPowVec_LFA(retStruct.instructright,:)),[.5 .5 1]);

        b(6)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.instructmouth,:),1),...
            SEM(zClusPowVec_LFA(retStruct.instructmouth,:)),[.7 .7 1]);

        
        %swagAxes(gca,fSize,'Time from INSTRUCT cue (ms)',['z power ' tStruct.retLbl1])
        swagAxes(gca,fSize,'','','Instruction')

        for bb = 1:length(b)
            set(b(bb),'facealpha',.75)
        end

          
        case 'moveWait' % this is move-wait 
        aLbl = 'mw';
        axes(thisAx.tf(3)); hold all    
    
            
        % HFA
        b(1)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.moveleft,:),1),...
            SEM(zClusPowVec_HFA(retStruct.moveleft,:)),[1 .3 .3]);
 
        b(2)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.moveright,:),1),...
            SEM(zClusPowVec_HFA(retStruct.moveright,:)),[1 .5 .5]);

        b(3)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_HFA(retStruct.movemouth,:),1),...
            SEM(zClusPowVec_HFA(retStruct.movemouth,:)),[1 .7 .7]);

        % LFA
        b(4)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.moveleft,:),1),...
            SEM(zClusPowVec_LFA(retStruct.moveleft,:)),[.3 .3 1]);
 
        b(5)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.moveright,:),1),...
            SEM(zClusPowVec_LFA(retStruct.moveright,:)),[.5 .5 1]);

        b(6)=plotSEM(tStruct.tBins,nanmean(zClusPowVec_LFA(retStruct.movemouth,:),1),...
            SEM(zClusPowVec_LFA(retStruct.movemouth,:)),[.7 .7 1]);
            
        
        %swagAxes(gca,fSize,'Time from GO cue (ms)',['z power ' tStruct.retLbl1])
        swagAxes(gca,fSize,'','','Move')
        
        %plot legend
        l = legend(b,{'HFA left','HFA right','HFA mouth','LFA left','LFA right','LFA mouth'},...
            'location','none','fontsize',10,'position',[.04 .85 0 0]);
        
         for bb = 1:length(b)
            set(b(bb),'facealpha',.75)
         end
         
    case 'moveInstruct'
        aLbl = 'mi';
            
    end
    
    %plot pow spectra
    freq_lim = [2 100];
    figure(h.(aLbl));
    thisAx.(aLbl)(5) = subplot(2,4,5); hold all
    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.instructleft),2)',...
    SEM(powSpectra(:,retStruct.instructleft)'),'r')

    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.waitleft),2)',...
    SEM(powSpectra(:,retStruct.waitleft)'),'g')

    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.moveleft),2)',...
    SEM(powSpectra(:,retStruct.moveleft)'),'b')
    swagAxes(gca,fSize,'Frequency (Hz)','Power','LEFT');set(gca,'xlim',freq_lim)
    legend('instruct','wait','move')

    thisAx.(aLbl)(6) = subplot(2,4,6); hold all    
    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.instructmouth),2)',...
    SEM(powSpectra(:,retStruct.instructmouth)'),'r')

    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.waitmouth),2)',...
    SEM(powSpectra(:,retStruct.waitmouth)'),'g')

    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.movemouth),2)',...
    SEM(powSpectra(:,retStruct.movemouth)'),'b')
    swagAxes(gca,fSize,'Frequency (Hz)','Power','MOUTH');set(gca,'xlim',freq_lim)
    legend('instruct','wait','move')
    
    thisAx.(aLbl)(7) = subplot(2,4,7); hold all
    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.instructright),2)',...
    SEM(powSpectra(:,retStruct.instructright)'),'r')

    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.waitright),2)',...
    SEM(powSpectra(:,retStruct.waitright)'),'g')

    plotSEM(config_pow_all.freQ(tStruct.fInd),...
    nanmean(powSpectra(:,retStruct.moveright),2)',...
    SEM(powSpectra(:,retStruct.moveright)'),'b')
    swagAxes(gca,fSize,'Frequency (Hz)','Power','RIGHT');set(gca,'xlim',freq_lim)
    legend('instruct','wait','move')
end


function [thisClusStruct] = getClusSelectivity(thisClusStruct,powStruct,config_pow_all,config_mn,plotFlag,thisAx,thisH,thisAx_tf,h_tf);


clusPow = shiftdim(nanmean(nanmean(powStruct.subjPow(thisClusStruct.freqIdx,thisClusStruct.timeIdx,:))),2);

% z-score across trials within this cluster
zClusPow = (clusPow - nanmean(clusPow))./nanstd(clusPow);

%% focus on frequency of interest for each category
zClusPowVec = shiftdim(nanmean(powStruct.subjPow(unique(thisClusStruct.freqIdx),:,:),1),1)'; 


%% Focus on time of interest
% power spectra during this cluster's time window of activation
% this will be useful when performing foof on these cluster activations
% later
clusPowSpectra = shiftdim(nanmean(powStruct.subjRawPow(:,thisClusStruct.timeIdx,:),2),2)';

% update clus pow struct
for f = 1:length(powStruct.retLbls)
    clusPowStruct.(powStruct.retLbls{f}) = zClusPow(powStruct.retStruct.(powStruct.retLbls{f}));
end

% %% implement  one-way t-stats as a fingerprint of cluster power
[~,ps ,~,stats] =cellfun(@ttest,struct2cell(clusPowStruct),'UniformOutput',false);
stats = [stats{:}];
thisClusStruct.ttestTstats = [stats.tstat];
thisClusStruct.ttestP = [ps{:}];

% means
meansClusPow =cellfun(@nanmean,struct2cell(clusPowStruct),'UniformOutput',false);
thisClusStruct.meansClusPow = [meansClusPow{:}];

% SEM
semsClusPow = cellfun(@SEM,struct2cell(clusPowStruct),'UniformOutput',false);
thisClusStruct.semsClusPow = [semsClusPow{:}];

% direction selectivity: one way anovas, t-test2, and zPowDiffVed
switch thisClusStruct.InstructOrMove
    case 0 % this is instruct-wait
    [p,tbl] = anovan([clusPowStruct.instructleft;clusPowStruct.instructright;clusPowStruct.instructmouth],...
    [1*ones(size(clusPowStruct.instructleft));2*ones(size(clusPowStruct.instructright));...
    3*ones(size(clusPowStruct.instructmouth))],'display','off');
    
    [~,p_rightleft,~,stats_rightleft] = ttest2(clusPowStruct.instructright,clusPowStruct.instructleft);
    [~,p_rightmouth,~,stats_rightmouth] = ttest2(clusPowStruct.instructright,clusPowStruct.instructmouth);
    [~,p_leftmouth,~,stats_leftmouth] = ttest2(clusPowStruct.instructleft,clusPowStruct.instructmouth);
    
    otherwise % this is a move-wait or move-instruct
    [p,tbl] = anovan([clusPowStruct.moveleft;clusPowStruct.moveright;clusPowStruct.movemouth],...
    [1*ones(size(clusPowStruct.moveleft));2*ones(size(clusPowStruct.moveright));...
    3*ones(size(clusPowStruct.movemouth))],'display','off');

    [~,p_rightleft,~,stats_rightleft] = ttest2(clusPowStruct.moveright,clusPowStruct.moveleft);
    [~,p_rightmouth,~,stats_rightmouth] = ttest2(clusPowStruct.moveright,clusPowStruct.movemouth);
    [~,p_leftmouth,~,stats_leftmouth] = ttest2(clusPowStruct.moveleft,clusPowStruct.movemouth);

end

% anova - direction selectivity
thisClusStruct.anova_dirSelectivity_p = p;
thisClusStruct.anova_dirSelectivity_F = tbl{2,6};
thisClusStruct.anova_dirSelectivity_df = tbl{2,3};
thisClusStruct.anova_dirSelectivity_MS = tbl{2,5};

% two way t-tests
thisClusStruct.ttest2_rightleft_t = stats_rightleft.tstat;
thisClusStruct.ttest2_rightleft_p = p_rightleft;

thisClusStruct.ttest2_rightmouth_t = stats_rightmouth.tstat;
thisClusStruct.ttest2_rightmouth_p = p_rightmouth;

thisClusStruct.ttest2_leftmouth_t = stats_leftmouth.tstat;
thisClusStruct.ttest2_leftmouth_p = p_leftmouth;

% update thisClusStruct
%thisClusStruct.clusPowStruct = clusPowStruct;
%thisClusStruct.zClusPowVec = zClusPowVec;


% local get Pow and plot function (wrapper)
function [tStruct_iw,tStruct_mw,tStruct_mi,config,h,thisAx] = getTstruct_local(subj,eLbl_list,config_mm, plotFlag, anatStruct)


[tStruct_iw,config] = le_calcTs_wrapper(subj,[],'motormap','instructWait',...
config_mm.fRange,config_mm.tRange_iw,config_mm.collapseFreqFlag,config_mm.powConfigNum,eLbl_list,config_mm.collapseElecFlag);

[tStruct_mw,config] = le_calcTs_wrapper(subj,[],'motormap','moveWait',...
config_mm.fRange,config_mm.tRange_mw,config_mm.collapseFreqFlag,config_mm.powConfigNum,eLbl_list,config_mm.collapseElecFlag);

[tStruct_mi,config] = le_calcTs_wrapper(subj,[],'motormap','moveInstruct',...
config_mm.fRange,config_mm.tRange_mw,config_mm.collapseFreqFlag,config_mm.powConfigNum,eLbl_list,config_mm.collapseElecFlag);

% update tStruct with anatomical struct
albls = fieldnames(anatStruct);
for f = 1:length(albls)
    switch albls{f}
        case {'anat','x','y','z','hemis','elecfullname','ROI'}
            tStruct_iw.(albls{f}) = anatStruct.(albls{f});
            tStruct_mw.(albls{f}) = anatStruct.(albls{f});
            tStruct_mi.(albls{f}) = anatStruct.(albls{f});
    end
end
tStruct_iw.uElbl = [tStruct_iw.subj '-' tStruct_iw.eLbl];
tStruct_mw.uElbl = [tStruct_mw.subj '-' tStruct_mw.eLbl];
tStruct_mi.uElbl = [tStruct_mi.subj '-' tStruct_mi.eLbl];



% 
thisAx = [];

if plotFlag
    h.tf = swagFig([-0.0789    0.0963    1.0789    0.7838]);
    set(gcf,'name',[anatStruct.eLbl '-' anatStruct.hemis '-' anatStruct.ROI '-TF'])
    for i = 1:8
        thisAx.tf(i) = subplot(2,4,i); 
    end    
    
    h.iw = swagFig([-0.0789    0.0963    1.0789    0.7838]);
    set(gcf,'name',[anatStruct.eLbl '-' anatStruct.hemis '-' anatStruct.ROI '-IW'])
    for i = 1:8
        thisAx.iw(i) = subplot(2,4,i); 
    end

    h.mw = swagFig([-0.0789    0.0963    1.0789    0.7838]);
    set(gcf,'name',[anatStruct.eLbl '-' anatStruct.hemis '-' anatStruct.ROI '-MW'])
    for i = 1:8
        thisAx.mw(i) = subplot(2,4,i); 
    end
    
    h.mi = swagFig([-0.0789    0.0963    1.0789    0.7838]);
    set(gcf,'name',[anatStruct.eLbl '-' anatStruct.hemis '-' anatStruct.ROI '-MI'])
    for i = 1:8
        thisAx.mi(i) = subplot(2,4,i); 
    end

    le_plotTF_local(tStruct_iw,tStruct_mw,tStruct_mi,config,0,0,thisAx);
else 
    h.tf = []; h.iw = [];h.mw = [];h.mi = [];
    thisAx.tf = [];thisAx.iw = []; thisAx.mw = [];thisAx.mi = [];
end



%% Config
function [ config_mm ] = getConfig_local(freqOption, compOption)
% This function is a config file w various parameters used throughout these
% analyses

if ~exist('freqOption','var') || isempty(freqOption)
    freqOption = 'all';
end
if ~exist('compOption','var') || isempty(compOption)
    compOption = 'moveWait';
end

switch compOption
    case 'moveWait'
        config_mm.comparison = 'moveWait';
 
    case 'instructWait'
        config_mm.comparison = 'instructWait';       
        
     case 'moveInstruct'
        config_mm.comparison = 'moveInstruct';          
end
%config_mm.comparison = 'instructWait';
%config_mm.comparison = 'rightWait';

% selectMotorSites
config_mm.hfa_range = [70 200];
config_mm.hfa_collapseFreqFlag = true;
config_mm.powConfigNum = 1;
config_mm.fdr_p_thresh = 0.05;
config_mm.num_sig_bins_thresh = 0; % this sets how many fdr-corrected sig. time windows you need to select a motor site

%freq def
config_mm.freqDef.HFA = [70 200];
config_mm.freqDef.LFA = [2 30];

% frequency bands of interest
config_mm.freqBandLims.HFA = [70 200];
config_mm.freqBandLims.gamma = [30 50];
config_mm.freqBandLims.beta = [12 30];
config_mm.freqBandLims.alpha = [8 12];
config_mm.freqBandLims.theta = [3 8];

% time bands of interest
config_mm.timeBandlims.prestim = [-500 0];
config_mm.timeBandlims.early = [0 500];
config_mm.timeBandlims.late = [500 1000];
config_mm.timeBandlims.sustained = [-500 1000];


% cluster motor sites based on HFA, tfPow and phase analyses
switch freqOption
    case 'all'
        config_mm.fRange = [3 200];
        config_mm.collapseFreqFlag = false;
    case 'hfa'
        config_mm.fRange = [70 200];
        config_mm.collapseFreqFlag = true;
    case 'gamma'
        config_mm.fRange = [30 50];
        config_mm.collapseFreqFlag = true;     
    case 'beta'
        config_mm.fRange = [12 30];
        config_mm.collapseFreqFlag = true; 
    case 'alpha'
        config_mm.fRange = [8 12];
        config_mm.collapseFreqFlag = true;   
    case 'theta'
        config_mm.fRange = [3 8];
        config_mm.collapseFreqFlag = true; 
end
config_mm.collapseElecFlag = false;

config_mm.collapseTimeFlag = true;
config_mm.tRange_mw = [-500 1000];
config_mm.tRange_iw = [-500 1000];
config_mm.tRange_mi = [-500 1000];


config_mm.tfplot_flag = 0; % sets whether to plot individual electrodes or not
config_mm.tfplot_printFlag = 1;

config_mm.phase_tRange = [-500 500];
config_mm.phase_collapseTimeFlag = 0;% [0 2000];

% clim for plot
config_mm.clim_tfDiff = [-7 7];

% plot preferences
config_mm.fSize = 30; 


function [] =  le_plotTF_local(tStruct_iw,tStruct_mw,tStruct_mi,config_pow,printFlag,printFolderLbl,thisAx);

%params
fSize = 22;

%make separate plots for Instruct-Wait (sensory) Move-Wait (motor)
%% Instruct/Wait
% get values for x-axis (time bins) and yaxis (freq bins, accounting for ...
%custom fRange used to generate this tStruct)
tBins = nanmean(config_pow.timeBins,2)'; 
tBins = tBins(tStruct_iw(1).tInd);
fBins = [config_pow.freqBins];
fBins = fBins(tStruct_iw(1).fInd); % filter by freq range used for this tStruct

% plot difference
% Z POW 1:  wait
axes(thisAx.iw(1)); hold all
imagesc(tBins,1:length(fBins),nanmean(tStruct_iw.pow2,3))
set(gca,'ydir','normal','clim',[-1 1])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'Time from GO cue (ms)','Frequency (Hz)',['zPow ' tStruct_iw.retLbl2])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')

% Z POW 2: Instruct
axes(thisAx.iw(2)); hold all
imagesc(tBins,1:length(fBins),nanmean(tStruct_iw.pow1,3))
set(gca,'ydir','normal','clim',[-1 1])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['zPow ' tStruct_iw.retLbl1])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')

% t(Z POW 1 - zPow 2): Instruct - Wait
axes(thisAx.tf(5)); hold all
meanT = tStruct_iw.tMat;
sigT = tStruct_iw.tMat_sig;
[~,thisH] = imagescnan(tBins,1:length(fBins),sigT);set(thisH,'visible','off');
thisH1 = imagescnan(tBins,1:length(fBins),meanT); alpha(.5)
set(gca,'ydir','normal','clim',[-7 7]);colormap('jet')
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['t(' tStruct_iw.retLbl1 '-' tStruct_iw.retLbl2 ')'])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
%c = colorbar('eastoutside');

% replicate in original figure
axes(thisAx.iw(3));hold all
copyobj(thisH1,thisAx.iw(3));alpha(1);
set(gca,'ydir','normal','clim',[-7 7])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['t(' tStruct_iw.retLbl1 '-' tStruct_iw.retLbl2 ')'])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')

%% plot MW
% get values for x-axis (time bins) and yaxis (freq bins, accounting for ...
%custom fRange used to generate this tStruct)
tBins = nanmean(config_pow.timeBins,2)'; 
tBins = tBins(tStruct_mw(1).tInd);
fBins = [config_pow.freqBins];
fBins = fBins(tStruct_mw(1).fInd); % filter by freq range used for this tStruct

% Z POW 1:  wait
axes(thisAx.mw(1)); hold all
imagesc(tBins,1:length(fBins),nanmean(tStruct_mw.pow2,3))
set(gca,'ydir','normal','clim',[-1 1])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'Time from GO cue (ms)','Frequency (Hz)',['zPow ' tStruct_mw.retLbl2])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')

% Z POW 1: move
axes(thisAx.mw(2)); hold all
imagesc(tBins,1:length(fBins),nanmean(tStruct_mw.pow1,3))
set(gca,'ydir','normal','clim',[-1 1])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['zPow ' tStruct_mw.retLbl1])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')

% t(Z POW 1 - zPow 2): Move - Wait
axes(thisAx.tf(6)); hold all
meanT = tStruct_mw.tMat;
sigT = tStruct_mw.tMat_sig;
[~,thisH] = imagescnan(tBins,1:length(fBins),sigT);set(thisH,'visible','off');
thisH1 = imagescnan(tBins,1:length(fBins),meanT); alpha(.5)
set(gca,'ydir','normal','clim',[-7 7]);colormap('jet')
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['t(' tStruct_mw.retLbl1 '-' tStruct_mw.retLbl2 ')'])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
%c = colorbar('eastoutside');
% Plot Y ticks for TF plots

% replicate in original figure
axes(thisAx.mw(3));hold all
copyobj(thisH1,thisAx.mw(3));alpha(1);
set(gca,'ydir','normal','clim',[-7 7])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['t(' tStruct_mw.retLbl1 '-' tStruct_mw.retLbl2 ')'])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')


%% plot MI (only difference
% get values for x-axis (time bins) and yaxis (freq bins, accounting for ...
%custom fRange used to generate this tStruct)
tBins = nanmean(config_pow.timeBins,2)'; 
tBins = tBins(tStruct_mi(1).tInd);
fBins = [config_pow.freqBins];
fBins = fBins(tStruct_mi(1).fInd); % filter by freq range used for this tStruct

% Z POW 1:  instruct
axes(thisAx.mi(1)); hold all
imagesc(tBins,1:length(fBins),nanmean(tStruct_mi.pow2,3))
set(gca,'ydir','normal','clim',[-1 1])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'Time from GO cue (ms)','Frequency (Hz)',['zPow ' tStruct_mw.retLbl2])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')

% Z POW 1: move
axes(thisAx.mi(2)); hold all
imagesc(tBins,1:length(fBins),nanmean(tStruct_mi.pow1,3))
set(gca,'ydir','normal','clim',[-1 1])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['zPow ' tStruct_mw.retLbl1])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')


% thresholding the t-statistics based on cluster stats
axes(thisAx.tf(7)); hold all
meanT = tStruct_mi.tMat;
sigT = tStruct_mi.tMat_sig;
[~,thisH] = imagescnan(tBins,1:length(fBins),sigT);set(thisH,'visible','off');
thisH1 = imagescnan(tBins,1:length(fBins),meanT); alpha(.5)
set(gca,'ydir','normal','clim',[-7 7]);colormap('jet')
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['t(' tStruct_mi.retLbl1 '-' tStruct_mi.retLbl2 ')'])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
%c = colorbar('location','manual','position',[.95 .55 0 0]);
% Plot Y ticks for TF plots

% replicate in original figure
axes(thisAx.mi(3));hold all
copyobj(thisH1,thisAx.mi(3));alpha(1);
set(gca,'ydir','normal','clim',[-7 7])
setFreqYTicks(gca,fBins)
swagAxes(gca,fSize,'','',['t(' tStruct_mi.retLbl1 '-' tStruct_mi.retLbl2 ')'])
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
c = colorbar('eastoutside');colormap('jet')


function [] = setFreqYTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'ytick',yt,...
'yticklabel',round(fBins(yt)))

function [] = updateTFplot_local(ax,thisClim,thisAnatStruct)

% plot colorbar
axes(ax(4)); 
set(gca,'clim',thisClim,'visible','off')
c = colorbar('location','westoutside');
set(c,'fontSize',16,'position',[.73 .1 .0158 .3]);

% update y lim on raw power plots
yL = [];
for i = 1:3
    axes(ax(i));
    yL = [yL get(gca,'ylim')];
end
for i = 1:3
    axes(ax(i));hold all
    set(gca,'ylim',[min(yL) max(yL)]);
    plot([0 0],get(gca,'ylim'),'linestyle','--','linewidth',3,'color', [.2 .2 .2])
end

% parse anatomical data
axes(ax(8))
set(gca,'position',[0.75    0.35    0.24    0.5]);hold all
isRightSide = [thisAnatStruct.x]>0;
if isRightSide
    plotBrain_local(1,[100 6]);
    %scatter3(thisAnatStruct.x,thisAnatStruct.y,thisAnatStruct.z,100,'filled','r')
    scatter3(100,thisAnatStruct.y,thisAnatStruct.z,100,'filled','r')

else
    plotBrain_local(1,[-100 6]);
    %scatter3(thisAnatStruct.x,thisAnatStruct.y,thisAnatStruct.z,100,'filled','r')
    scatter3(-100,thisAnatStruct.y,thisAnatStruct.z,100,'filled','r')

end


function [h] = plotBrain_local(facealpha,viewAngle)

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
camzoom(1.5)
% GRAVEYARD %


% % % aggregate clusStruct - % only use if multiple cluster based comparisons
% if ~exist('clusStruct_iw','var') || isempty(fieldnames(clusStruct_iw))
%     clusStruct = clusStruct_mw;
% elseif isempty(fieldnames(clusStruct_mw))
%     clusStruct = clusStruct_iw;
% else
%     clusStruct = [clusStruct_iw clusStruct_mw];
% end
% 
% % sort the activations by largest cluster statistic
% if ~isempty(fieldnames(clusStruct))
%     [~,sortIdx] = sort([clusStruct.stat]);
%     clusStruct(sortIdx);
% end

