function [statsStruct] = le_mn_ROI(plotFlag,printFlag)
% This function performs an ROI analysis. Performs the following:
% 1) collects all significant activations by subj
% 2) plots sig activations on brain 
% 2) categorizes activations into high/low freq; pos/neg for whole brain and ROI
% 3) freq characteristics; increases and decreases; whole brain v roi
% 4) do tilt type effects co-occur at the same site? whole brain v roi
% 5) timing analysis
% 6) movt-selectivity; whole brain v roi

if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = 1;
end
if ~exist('printFlag','var') || isempty(printFlag)
    printFlag = 0;
end
%% initialize
dirs = le_dirs;

params.subj_list = {'HUP084','HUP087','HUP088','HUP089','HUP090','HUP090_1',...
    'HUP091','HUP092','HUP111','HUP112'};
params.saveDir = fullfile(dirs.scratch,'clusStruct');
params.saveFile = ['motornetwork-groupClusStruct-' params.subj_list{:} '.mat'];
params.printDir = fullfile(dirs.scratch,'figs','le_mn_ROI',date);
params.statsDir = fullfile(dirs.scratch,'stats','le_mn');

% print
params.fileformat = '-depsc2';
params.resolution = '-r300';%-r0

% cluster effect filtering
% anatomy
% anatomy
params.filt_roi = {'all','perirolandic','prefrontal','temporal','parietal','MTL','occipital'};
params.filt_hemis = {'','','','','','',''};
params.filt_subj = {''};
params.filt_roi_kol= {'r','g'};%;'b','y','m','c'};%

for i =1:length(params.filt_roi)
    params.filt_roiLbl{i}= [params.filt_hemis{i} params.filt_roi{i}];
end

% % freq band 
params.freqBandLims.HFA = [70 200];
params.freqBandLims.gamma = [30 50];
params.freqBandLims.beta = [12 30];
params.freqBandLims.alpha = [10 12];
params.freqBandLims.theta = [3 8];

% time bands of interest
params.timeBandlims.prestim = [-500 0];
params.timeBandlims.early = [0 500];
params.timeBandlims.late = [500 1000];
params.timeBandlims.sustained = [-500 1000];

% selectivity
params.pThresh = 0.01;
params.pThresh_moveSelective = 0.01;

%clus stats
params.clusConfig.alfa = [0.025 0.975];
params.clusConfig.iters =10;%1000

% fig params
% TF plot
params.fig.tf.fSize = 12;
params.fig.tf.subplot_x = 1;
params.fig.tf.subplot_y = 3;
params.fig.tf.clim = [-7 7];
params.fig.tf.collapseWinSubj = 0;

% zPow plot
params.fig.zpow.fSize = 12;
params.fig.zpow.subplot_x = 1;
params.fig.zpow.subplot_y = 3;
params.fig.zpow.ylim = [-.5 .32];
params.fig.zpow.collapseWinSubj = 0;

% zCorr plot
params.fig.zcorr.fSize = 12;
params.fig.zcorr.subplot_x = 1;
params.fig.zcorr.subplot_y = 2;

% counts
params.fig.counts.fSize = 12;
params.fig.counts.toiLblsToInclude = {'prestim','early','late','sustained'};
params.fig.counts.foiLblsToInclude = {'HFA','beta','theta'};
params.fig.counts.doSubplot = 1;

params.fig.counts.subplot_x = length(params.fig.counts.foiLblsToInclude);;
params.fig.counts.subplot_y = 3;

%% make example figures

% instruction-related HFA example
[tStruct_iw_eg1] = le_mn_spectralChange('HUP087','20-21',0,0);

% movement-related HFA example
[~,tStruct_mw_eg2] = le_mn_spectralChange('HUP087','22-23',0,0);

% collect p's and save in stats dir

% load pStruct
pStruct = loadPStruct_local(params);

% overwrite example fields
pStruct.eg_instruct = [tStruct_iw_eg1.posClusStatsP tStruct_iw_eg1.negClusStatsP];
pStruct.eg_move = [tStruct_mw_eg2.posClusStatsP tStruct_mw_eg2.negClusStatsP];

% save updated pStruct and clear
save('pStruct','pStruct');  clear pStruct;

%% load clusStruct group for all subjects
[tStruct_iw_grp,tStruct_mw_grp,tStruct_mi_grp,...
    clusStruct_iw_grp,clusStruct_mw_grp,clusStruct_mi_grp] = ...
    loadGroupClusStruct_local(params.subj_list,params.saveDir,params.saveFile);
disp('group data LOADED...Ready to Go!!!')

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


%% TF plot (grand average of sig electrodes within a region)
[h,statsStruct] = tfPlot_local(roiStruct,params);

%% count cluster effects (by any frequency; any time)
params.fig.counts.toiLblsToInclude = {'all'};
params.fig.counts.foiLblsToInclude = {'all','HFA','gamma','beta','alpha','theta'};
params.fig.counts.doSubplot = 0;
[hh, statsStruct, countStruct, expStruct] = doCounts_local(roiStruct,params,statsStruct);
h = [h hh];

%% stats counts
[statsStruct,structForTable] = doStats_counts_local(countStruct,expStruct,params,statsStruct);

% do FDR correction
%fdr(pVec,0.001);


%% count cluster effects (by frequency; any time)
params.fig.counts.toiLblsToInclude = {'all'};
params.fig.counts.foiLblsToInclude = {'HFA','gamma','beta','alpha','theta'};
params.fig.counts.doSubplot = 0;
[hh, statsStruct] = doCounts_local(roiStruct,params,statsStruct);
h = [h hh];


%% pow comparison (within region)
[hh,statsStruct] = zPowPlot_local(roiStruct,params,statsStruct);
h = [h hh];
%% pow correlations (across electrodes)
[hh,statsStruct] = zPowCorr_local(roiStruct,params,statsStruct);
h = [h hh];
%% count cluster effects (by time of intrest) (inc vs. decreases, effector selectivity)
params.fig.counts.toiLblsToInclude = {'prestim','early','late','sustained'};
params.fig.counts.foiLblsToInclude = {'HFA','beta','theta'};
params.fig.counts.doSubplot =  1;
[hh, statsStruct] = doCounts_local(roiStruct,params,statsStruct);
h = [h hh];
%% count cluster effects (by frequency; any time)
params.fig.counts.toiLblsToInclude = {'all'};
params.fig.counts.foiLblsToInclude = {'HFA','gamma','beta','alpha','theta'};
params.fig.counts.doSubplot = 0;
[hh, statsStruct] = doCounts_local(roiStruct,params,statsStruct);
h = [h hh];

%% do stats
%keyboard
%% print
if printFlag
    cd_mkdir(params.printDir)
    for i = 1:length(h)
        figure(h(i));
        print(h(i),get(h(i),'Name'),params.fileformat,params.resolution);
        pause(3)
        close
    end
end

%%
function[statsStruct,structForTable] = doStats_counts_local(countStruct,expStruct,params,statsStruct)


%% Power Changes v.s. Chance, whole brain

% toi is the time of interest (default = 'all')
toiLbls =  params.fig.counts.toiLblsToInclude;

for t = 1:length(toiLbls)
toi = toiLbls{t};%

    % foi is the frequency label of interest
    foiLbls = params.fig.counts.foiLblsToInclude;
    for f = 1:length(foiLbls)
        foi = foiLbls{f};
        for i = 1:length(params.filt_roi)
            % choose region
            roi = params.filt_roi{i};% 'all';

            %Instruction vs. Wait

            % initialize counts
            %all regions, observed sign. increases and decreases, and expected by
            %chance (set by pThreshold)
            totCount = countStruct.(roi).IW.(foi).(toi).total;%total count
            sigIncrease_o = countStruct.(roi).IW.(foi).(toi).sigIncrease;
            sigDecrease_o = countStruct.(roi).IW.(foi).(toi).sigDecrease;
            effSel_o = countStruct.(roi).IW.(foi).(toi).effectorSelective;
            sigIncrease_e = expStruct.(roi).IW.(foi).(toi).sigIncrease;
            sigDecrease_e = expStruct.(roi).IW.(foi).(toi).sigDecrease;
            effSel_e = expStruct.(roi).IW.(foi).(toi).effectorSelective;

            % do we observe more significant increases as compared to chance?
            % expected frequencies determined by false positive rate of
            % pThresh (e.g., 0.01) per electrode
            [statsStruct.countsIncreaseIW.(roi).chi2,statsStruct.countsIncreaseIW.(roi).p,...
                statsStruct.countsIncreaseIW.(roi).df] = ...
                chi2test([sigIncrease_o totCount-sigIncrease_o;sigIncrease_e totCount-sigIncrease_e]);

            % do we observe more significant decreases as compared to chance?
            [statsStruct.countsDecreaseIW.(roi).chi2,statsStruct.countsDecreaseIW.(roi).p,...
                statsStruct.countsDecreaseIW.(roi).df] = ...
                chi2test([sigDecrease_o totCount-sigDecrease_o;sigDecrease_e totCount-sigDecrease_e]);

            % Unclear what the expected frequency of effector selectivity
            % per electrode would be. Leave this out for now
%             %do we observe more effector selectivity as compared to chance?
%             [statsStruct.effSelIW.(roi).chi2,statsStruct.effSelIW.(roi).p,...
%             statsStruct.effSelIW.(roi).df] = ...
%             chi2test([effSel_o totCount-effSel_o;effSel_e totCount-effSel_e]);
% 

            %Move vs. Wait
            % initialize counts
            %all regions, observed sign. increases and decreases, and expected by
            %chance (set by pThreshold)
            totCount = countStruct.(roi).MW.(foi).(toi).total;%total count
            sigIncrease_o = countStruct.(roi).MW.(foi).(toi).sigIncrease;
            sigDecrease_o = countStruct.(roi).MW.(foi).(toi).sigDecrease;
            effSel_o = countStruct.(roi).MW.(foi).(toi).effectorSelective;
            sigIncrease_e = expStruct.(roi).MW.(foi).(toi).sigIncrease;
            sigDecrease_e = expStruct.(roi).MW.(foi).(toi).sigDecrease;
            effSel_e = expStruct.(roi).MW.(foi).(toi).effectorSelective;

            % do we observe more significant increases as compared to chance?
            [statsStruct.countsIncreaseMW.(roi).chi2,statsStruct.countsIncreaseMW.(roi).p,...
                statsStruct.countsIncreaseMW.(roi).df] = ...
                chi2test([sigIncrease_o totCount-sigIncrease_o;sigIncrease_e totCount-sigIncrease_e]);

            % do we observe more significant decreases as compared to chance?
            [statsStruct.countsDecreaseMW.(roi).chi2,statsStruct.countsDecreaseMW.(roi).p,...
                statsStruct.countsDecreaseMW.(roi).df] = ...
                chi2test([sigDecrease_o totCount-sigDecrease_o;sigDecrease_e totCount-sigDecrease_e]);

            % see above note
%            % more effector selectivitvity than chance?
%             [statsStruct.effSelMW.(roi).chi2,statsStruct.effSelMW.(roi).p,...
%             statsStruct.effSelMW.(roi).df] = ...
%             chi2test([effSel_o totCount-effSel_o;effSel_e totCount-effSel_e]);
% 

           % make struct for table
           structForTable.(foi).(toi).counts_IW(i) = countStruct.(roi).IW.(foi).(toi);
           structForTable.(foi).(toi).counts_MW(i) = countStruct.(roi).MW.(foi).(toi);
           structForTable.(foi).(toi).chi2Increase_IW(i)= statsStruct.countsIncreaseIW.(roi);
           structForTable.(foi).(toi).chi2Decrease_IW(i)= statsStruct.countsDecreaseIW.(roi);
           structForTable.(foi).(toi).chi2Increase_MW(i)= statsStruct.countsIncreaseMW.(roi);
           structForTable.(foi).(toi).chi2Decrease_MW(i)= statsStruct.countsDecreaseMW.(roi);
           %structForTable.(foi).(toi).chi2EffSel_IW(i)= statsStruct.effSelIW.(roi);
           %structForTable.(foi).(toi).chi2EffSel_MW(i)= statsStruct.effSelMW.(roi);

        end
    end
end

% Collect p-values for chi-2 tests comparing frequency of increases and
% decreases vs. chance (no freq of interest or  or time of selection)
foiLbl = 'all';
toiLbl = 'all';
tblLbl_list = fieldnames(structForTable.(foiLbl).(toiLbl));

% load pStruct (will overwrite this w p-values from chi 2 test)
% will be used later for fdr correction in le_mn_fdr
pStruct = loadPStruct_local(params);

% pVec to hold these p-values, will update pStruct with these values below
% FDR correction
pVec = [];
for i = 1:length(tblLbl_list)
    tblLbl = tblLbl_list{i};
    T.(tblLbl) = struct2table(structForTable.(foiLbl).(toiLbl).(tblLbl),'RowNames',params.filt_roi);
    
    % populate pVec
    switch tblLbl
        case{'chi2Increase_IW','chi2Decrease_IW','chicl2Increase_MW','chi2Decrease_MW',...
                'chi2EffSel_IW','chi2EffSel_MW'}
            pVec = [pVec [structForTable.(foiLbl).(toiLbl).(tblLbl).p]];
    end
    
%     if printFlag == 1
%         cd_mkdir(params.printDir)
%         writetable(T.(tblLbl),fullfile(params.printDir,tblLbl),'writerownames',true,'delimiter',' ');
%     end
end

% update pStruct
pStruct.powChangesVsChance = pVec;
save('pStruct','pStruct'); clear pStruct;



%% Next, we compute the fruency of power changes per frequency band.
toiLbl = 'all';
% do frequency-based counts
incMatIW_freq = nan(length(params.filt_roi), length(params.fig.counts.foiLblsToInclude));
decMatIW_freq = incMatIW_freq;
incMatMW_freq = incMatIW_freq;
decMatMW_freq = incMatIW_freq;
effSelMatIW_freq = incMatIW_freq;
effSelMatMW_freq = incMatIW_freq;
totCounts_freq = incMatIW_freq;

for i = 1:length(params.fig.counts.foiLblsToInclude)
    foiLbl = params.fig.counts.foiLblsToInclude{i};
    
    totCounts_freq(:,i) = [structForTable.(foiLbl).(toiLbl).counts_IW.total]';
    incMatIW_freq(:,i) = [structForTable.(foiLbl).(toiLbl).counts_IW.sigIncrease]';
    decMatIW_freq(:,i) = [structForTable.(foiLbl).(toiLbl).counts_IW.sigDecrease]';
    effSelMatIW_freq(:,i) = [structForTable.(foiLbl).(toiLbl).counts_IW.effectorSelective]';

    
    incMatMW_freq(:,i) = [structForTable.(foiLbl).(toiLbl).counts_MW.sigIncrease]';
    decMatMW_freq(:,i) = [structForTable.(foiLbl).(toiLbl).counts_MW.sigDecrease]';
    effSelMatMW_freq(:,i) = [structForTable.(foiLbl).(toiLbl).counts_MW.effectorSelective]';
end

% convert to percentages
incMatIW_freq_prct =100*(incMatIW_freq./totCounts_freq);
decMatIW_freq_prct =100*(decMatIW_freq./totCounts_freq);
effSelMatIW_freq_prct =100*(effSelMatIW_freq./totCounts_freq);


incMatMW_freq_prct =100*(incMatMW_freq./totCounts_freq);
decMatMW_freq_prct =100*(decMatMW_freq./totCounts_freq);
effSelMatMW_freq_prct =100*(effSelMatMW_freq./totCounts_freq);


%% additional stats for freq band, regional and move-sel effects:

pStruct = loadPStruct_local(params);
pVec = [];

%total number of movement-specific changes in instruct vs. move
iw_effSel_all = structForTable.all.all.counts_IW(1).effectorSelective;
mw_effSel_all = structForTable.all.all.counts_MW(1).effectorSelective;
total =  structForTable.all.all.counts_IW(1).total;

[chi2,p,df] = chi2test([iw_effSel_all total-iw_effSel_all;mw_effSel_all total-total]);
pVec = [pVec p];

% direct comparison between regional distribution 
% of power increases instruct vs. move
a = T.counts_IW{2:end,2}; % this measures the distribution of instruction related INCREASES across regions
b = T.counts_MW{2:end,2}; % this measures the distribution of movement related INCREASES across regions
[chi2,p,df] = chi2test([a';b']);
pVec = [pVec p];

% direct comparison between regional distribution 
% of power decreases instruct vs. move
a = T.counts_IW{2:end,4}; % this measures the distribution of instruction related DECREASES across regions
b = T.counts_MW{2:end,4}; % this measures the distribution of movement related DECREASES across regions
[chi2,p,df] = chi2test([a';b']);
pVec = [pVec p];

% do HFA, beta and theta show differences in the direction of power change
% across the brain?

all_changes = [countStruct.all.IW.all.all.sigIncrease,...
    countStruct.all.IW.all.all.sigDecrease,...
    countStruct.all.MW.all.all.sigIncrease,...
    countStruct.all.MW.all.all.sigDecrease];

HFA_changes = [countStruct.all.IW.HFA.all.sigIncrease,...
    countStruct.all.IW.HFA.all.sigDecrease,...
    countStruct.all.MW.HFA.all.sigIncrease,...
    countStruct.all.MW.HFA.all.sigDecrease]; % theta dist of power changes

gamma_changes = [countStruct.all.IW.gamma.all.sigIncrease,...
    countStruct.all.IW.gamma.all.sigDecrease,...
    countStruct.all.MW.gamma.all.sigIncrease,...
    countStruct.all.MW.gamma.all.sigDecrease]; % theta dist of power changes


beta_changes = [countStruct.all.IW.beta.all.sigIncrease,...
    countStruct.all.IW.beta.all.sigDecrease,...
    countStruct.all.MW.beta.all.sigIncrease,...
    countStruct.all.MW.beta.all.sigDecrease]; % theta dist of power changes

alpha_changes = [countStruct.all.IW.alpha.all.sigIncrease,...
    countStruct.all.IW.alpha.all.sigDecrease,...
    countStruct.all.MW.alpha.all.sigIncrease,...
    countStruct.all.MW.alpha.all.sigDecrease]; % theta dist of power changes

theta_changes = [countStruct.all.IW.theta.all.sigIncrease,...
    countStruct.all.IW.theta.all.sigDecrease,...
    countStruct.all.MW.theta.all.sigIncrease,...
    countStruct.all.MW.theta.all.sigDecrease]; % theta dist of power changes
[chi2,p,df] = chi2test([HFA_changes;all_changes]);
pVec = [pVec p];
[chi2,p,df] = chi2test([gamma_changes;all_changes]);
pVec = [pVec p];
[chi2,p,df] = chi2test([beta_changes;all_changes]);
pVec = [pVec p];
[chi2,p,df] = chi2test([alpha_changes;all_changes]);
pVec = [pVec p];
[chi2,p,df] = chi2test([theta_changes;all_changes]);
pVec = [pVec p];
[chi2,p,df] = chi2test([theta_changes;beta_changes]);
pVec = [pVec p];

% save pStruct
save('pStruct','pStruct'); clear pStruct

% % (   ) Currently not inculded: 
% % The below code compares the distribution of various frequency band
% % changes across regions to assess whether there are regional differences
% % in frequency bands. 
% tot_counts = totCounts_freq(2:end,1); % all electrodes in each of the regions
% %total HFA counts in all regions
% HFA_counts = incMatIW_freq(2:end,2)+decMatIW_freq(2:end,2)+...
%     incMatMW_freq(2:end,2)+decMatMW_freq(2:end,2);
% %total gamma counts in all regions
% gamma_counts = incMatIW_freq(2:end,3)+decMatIW_freq(2:end,3)+...
%     incMatMW_freq(2:end,3)+decMatMW_freq(2:end,3);
% %total beta counts in all regions
% beta_counts = incMatIW_freq(2:end,4)+decMatIW_freq(2:end,4)+...
%     incMatMW_freq(2:end,4)+decMatMW_freq(2:end,4);
% %total alpha counts in all regions
% alpha_counts = incMatIW_freq(2:end,5)+decMatIW_freq(2:end,5)+...
%     incMatMW_freq(2:end,5)+decMatMW_freq(2:end,5);
% %total theta counts in all regions
% theta_counts = incMatIW_freq(2:end,6)+decMatIW_freq(2:end,6)+...
%     incMatMW_freq(2:end,6)+decMatMW_freq(2:end,6);
% 
% % HFA vs. total
% [chi2,p,df] = chi2test([HFA_counts';tot_counts']);
% pVec = [pVec p];
% % gamma vs. total
% [chi2,p,df] = chi2test([gamma_counts';tot_counts']);
% pVec = [pVec p];
% % beta vs. total
% [chi2,p,df] = chi2test([beta_counts';tot_counts']);
% pVec = [pVec p];
% % alpha vs. total
% [chi2,p,df] = chi2test([alpha_counts';tot_counts']);
% pVec = [pVec p];
% % theta vs. total
% [chi2,p,df] = chi2test([theta_counts';tot_counts']);
% pVec = [pVec p];
% 

function [h, statsStruct,countStruct,expStruct] = doCounts_local(specStruct_ROI,params,statsStruct);
% loop through regions
for i = 1:length(params.filt_roi)
    % get vars for this ROI
    roi = params.filt_roiLbl{i};
    roiStruct = specStruct_ROI.(roi);
    % create figure for each region
    if params.fig.counts.doSubplot
        h(i) = swagFig([ 0    0.1512    0.9602    0.7288]);
        h(i).Name = ['zCountsByTime-' roi];
    else
        h(i) = swagFig([0    0.4975    1.0000    0.3825]);
        h(i).Name = ['zCounts-' roi];
    end

    
    % make subplots
    % each region has a figure. Each freq has a subplot (row)
    % each subplot has a set of bar graphs for each time of interst
    % (prestim, early, late, sustained)
    % first bar shows % sig electrodes with activation in that FOI/TOI (%
    % retained)
    % second bar shows % of FOI/TOI sig electrodes with increases
    % third bar shows % of FOI/TOI sig electrodes with decreases
    % fourth bar shows % of FOI/TOI sig electrodes with effector
    % selectivity
    % ~~~~ probably need bars of different types w legend
    subplot_i = 0;
    foiLbls = params.fig.counts.foiLblsToInclude;
    toiLbls = params.fig.counts.toiLblsToInclude;
    countMat_prct_iw = cell(1,length(foiLbls));% nan(length(foiLbls),length(toiLbls));
    countMat_prct_mw = cell(1,length(foiLbls));% nan(length(foiLbls),length(toiLbls));
    countMat_prct_mi = cell(1,length(foiLbls));%nan(length(foiLbls),length(toiLbls));
    for j = 1:length(foiLbls)
        % Instruct - wait 
        subplot_i = subplot_i+1;
        compLbl = 'instructWait';
        foiLbl =foiLbls{j}; % freq of interest label
        [countStruct.(roi).IW.(foiLbl),countMat_prct_iw{j},expStruct.(roi).IW.(foiLbl)] = doCountsFOI_local(roiStruct,compLbl,foiLbl,params,roi,subplot_i);

        % MOVE-Wait 
        subplot_i = subplot_i+1;
        i = i+1;
        compLbl = 'moveWait';
        [countStruct.(roi).MW.(foiLbl),countMat_prct_mw{j},expStruct.(roi).MW.(foiLbl)] = doCountsFOI_local(roiStruct,compLbl,foiLbl,params,roi,subplot_i);

        %  MOVE - Instruct
        subplot_i = subplot_i+1;
        compLbl = 'moveInstruct';
        [countStruct.(roi).MI.(foiLbl),countMat_prct_mi{j},expStruct.(roi).MI.(foiLbl)] = doCountsFOI_local(roiStruct,compLbl,foiLbl,params,roi,subplot_i);

        %% stats
        % START HERE
        %keyboard
        % expected counts vs actual counts
        
        
       
        
    end
    if length(foiLbls)==1&&length(toiLbls)==1
        subplot(1,3,1)
        countMat_prct_iw_mat = cat(1,countMat_prct_iw{:});
        hold all
        b = bar(countMat_prct_iw_mat,'linewidth',1.5);
        swagAxes(gca,params.fig.counts.fSize,'',...
            '% all electrodes in region',[roi ' instruct - wait (n electrodes = ' num2str(length(roiStruct.tStruct_iw_grp)) ')']);
        set(gca,'ylim',[0 50],'xlim',[.5 3.5],'xtick',1:3,'xticklabel',{'Increase','Decrease','Movement Selective'})
      
        subplot(1,3,2)
        countMat_prct_mw_mat = cat(1,countMat_prct_mw{:});
        hold all
        b = bar(countMat_prct_mw_mat,'linewidth',1.5);
        swagAxes(gca,params.fig.counts.fSize,'',...
            '% all electrodes in region',[roi ' move - wait (n electrodes = ' num2str(length(roiStruct.tStruct_mw_grp)) ')']);
        set(gca,'ylim',[0 50],'xlim',[.5 3.5],'xtick',1:3,'xticklabel',{'Increase','Decrease','Movement Selective'})
  
        subplot(1,3,3)
        countMat_prct_mi_mat = cat(1,countMat_prct_mi{:});
        hold all
        b = bar(countMat_prct_mi_mat,'linewidth',1.5);
        swagAxes(gca,params.fig.counts.fSize,'',...
            '% all electrodes in region',[roi ' move - instruct (n electrodes = ' num2str(length(roiStruct.tStruct_mi_grp)) ')']);
        set(gca,'ylim',[0 50],'xlim',[.5 3.5],'xtick',1:3,'xticklabel',{'Increase','Decrease','Effector Selective'})
    elseif params.fig.counts.doSubplot == 0
        subplot(1,3,1)
        countMat_prct_iw_mat = cat(1,countMat_prct_iw{:});
        hold all
        b = bar(countMat_prct_iw_mat,'linewidth',1.5);
        swagAxes(gca,params.fig.counts.fSize,'',...
            '% all electrodes in region',[roi ' instruct - wait (n electrodes = ' num2str(length(roiStruct.tStruct_iw_grp)) ')']);
        set(gca,'ylim',[0 50],'xlim',[.5 length(foiLbls)+.5],'xtick',1:length(foiLbls),'xticklabel',foiLbls)
      
        subplot(1,3,2)
        countMat_prct_mw_mat = cat(1,countMat_prct_mw{:});
        hold all
        b = bar(countMat_prct_mw_mat,'linewidth',1.5);
        swagAxes(gca,params.fig.counts.fSize,'',...
            '% all electrodes in region',[roi ' move - wait (n electrodes = ' num2str(length(roiStruct.tStruct_mw_grp)) ')']);
        set(gca,'ylim',[0 20],'xlim',[.5 length(foiLbls)+.5],'xtick',1:length(foiLbls),'xticklabel',foiLbls)
  
        subplot(1,3,3)
        countMat_prct_mi_mat = cat(1,countMat_prct_mi{:});
        hold all
        b = bar(countMat_prct_mi_mat,'linewidth',1.5);
        swagAxes(gca,params.fig.counts.fSize,'',...
            '% all electrodes in region',[roi ' move - instruct (n electrodes = ' num2str(length(roiStruct.tStruct_mi_grp)) ')']);
        set(gca,'ylim',[0 20],'xlim',[.5 length(foiLbls)+.5],'xtick',1:length(foiLbls),'xticklabel',foiLbls)
       % l = legend(b,fieldnames(countStruct.MI_HFA.(toiLbls{1})),'location','northeast');%,'location',[0.85 0.2093 0.0720 0.0986]);
    end  
    
    
    
end

function [countStruct,countMat_prct,expStruct] = doCountsFOI_local(roiStruct,compLbl,foiLbl,params,roi,subplot_i)
switch compLbl 
    case 'instructWait'
        thisClusStruct = roiStruct.clusStruct_iw_grp;
        numElec_tot = length(roiStruct.tStruct_iw_grp);
    case 'moveWait'
        thisClusStruct = roiStruct.clusStruct_mw_grp;
        numElec_tot = length(roiStruct.tStruct_mw_grp);
    case 'moveInstruct'
        thisClusStruct = roiStruct.clusStruct_mi_grp;
        numElec_tot = length(roiStruct.tStruct_mi_grp);
end

toiLbls = params.fig.counts.toiLblsToInclude;
for t = 1:length(toiLbls)
    % prestim v early v late v sustained v any
    [countStruct.(toiLbls{t}),expStruct.(toiLbls{t})]...
        = doCountsFOITOI_local(foiLbl,toiLbls{t},thisClusStruct,params,numElec_tot);
end


%this sfcn will take a structure and plot a group of bars)
flbls = fieldnames(countStruct);
fflbls = fieldnames(countStruct.(flbls{1}));
countMat = nan(length(flbls),length(fflbls));
for f = 1:length(flbls)
 c = struct2cell(countStruct.(flbls{f}));
 countMat(f,:) = [c{:}];
end

countMat_prct = (countMat./numElec_tot)*100;

if params.fig.counts.doSubplot
    subplot(params.fig.counts.subplot_x,params.fig.counts.subplot_y,subplot_i)
    hold all
    b = bar(countMat_prct,'linewidth',1.5);
    swagAxes(gca,params.fig.counts.fSize,'',...
        '% all electrodes in region',[roi ' ' foiLbl ' ' compLbl ' (n electrodes = ' num2str(numElec_tot) ')']);
    set(gca,'ylim',[0 25],'xlim',[.5 length(flbls)+.5],'xtick',1:length(flbls),'xticklabel',flbls)
    
    if subplot_i == 1
        l = legend(b,fflbls,'location','northwest');%,'location',[0.85 0.2093 0.0720 0.0986]);
    end

end
%legend(b,fflbls,'Location','northwest')
function[thisCountsStruct, thisExpStruct]...
    = doCountsFOITOI_local(foi,toi,thisClusStruct,params,numElec_tot);
%thisCountsStruct.numSigElec_tot = length(unique({thisClusStruct.uElbl}));

if strcmp(toi,'all') % if all time bins are selected, it will count based on foi alone
    if strcmp(foi,'all')
        retIdx = true(size(thisClusStruct));
    else
        retIdx = strcmp(foi,{thisClusStruct.freqBandLbl});
    end
else
    retIdx = strcmp(foi,{thisClusStruct.freqBandLbl})&...
        strcmp(toi,{thisClusStruct.timeBandLbl});
end
%thisCountsStruct.significant = length(unique({thisClusStruct(retIdx).uElbl}));
thisCountsStruct.total = numElec_tot;
thisCountsStruct.sigIncrease = sum([thisClusStruct(retIdx).clusPolarity]==1&[thisClusStruct(retIdx).p]<=params.pThresh);
thisCountsStruct.sigIncrease_prct = (thisCountsStruct.sigIncrease/numElec_tot)*100;
thisCountsStruct.sigDecrease = sum([thisClusStruct(retIdx).clusPolarity]==-1&[thisClusStruct(retIdx).p]<=params.pThresh);
thisCountsStruct.sigDecrease_prct = (thisCountsStruct.sigDecrease/numElec_tot)*100;
thisCountsStruct.effectorSelective= sum([thisClusStruct(retIdx).anova_dirSelectivity_p]<=params.pThresh_moveSelective);
thisCountsStruct.effectorSelective_prct = (thisCountsStruct.effectorSelective/numElec_tot)*100;

%thisCountsStruct.effectorSelectiveIncrease = sum([thisClusStruct(retIdx).anova_dirSelectivity_p]<=params.pThresh_moveSelective&[thisClusStruct(retIdx).clusPolarity]==1);
%thisCountsStruct.effectorSelectiveDecrease = sum([thisClusStruct(retIdx).anova_dirSelectivity_p]<=params.pThresh_moveSelective&[thisClusStruct(retIdx).clusPolarity]==-1);

thisExpStruct.sigIncrease = numElec_tot*params.pThresh;
thisExpStruct.sigDecrease = numElec_tot*params.pThresh;
thisExpStruct.effectorSelective = numElec_tot*params.pThresh_moveSelective;

function [h,statsStruct] = zPowCorr_local(specStruct_ROI,params,statsStruct)
tBins = specStruct_ROI.all.tStruct_iw_grp.tBins;
fBins = specStruct_ROI.all.tStruct_iw_grp.fBins;
fSize = params.fig.zcorr.fSize;

% get time bins of interest
[~,tStart_early] = min(abs(tBins - params.timeBandlims.early(1)));
[~,tEnd_early] = min(abs(tBins - params.timeBandlims.early(2)));

[~,tStart_late] = min(abs(tBins - params.timeBandlims.late(1)));
[~,tEnd_late] = min(abs(tBins - params.timeBandlims.late(2)));


% loop through regions
for i = 1:length(params.filt_roi)
    % get vars for this ROI
    roi = params.filt_roiLbl{i};

    % create figure for each region
    h(i) = swagFig([0    0.5350    0.5773    0.3450]);
    h(i).Name = ['zPowCorr-' roi];
    

    % skip this if no activations
    if ~isfield(specStruct_ROI,params.filt_roiLbl{i})
        continue
    end

    % collect z pow variables
    thisHFA_wait = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zHFA_wait_mean);
    thisHFA_instruct = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zHFA_instruct_mean);
    thisHFA_move = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zHFA_move_mean);

    thisBeta_wait = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zBeta_wait_mean);
    thisBeta_instruct = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zBeta_instruct_mean);
    thisBeta_move = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zBeta_move_mean);

    thisTheta_wait = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zTheta_wait_mean);
    thisTheta_instruct = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zTheta_instruct_mean);
    thisTheta_move = cat(1,specStruct_ROI.(roi).tStruct_iw_grp.zTheta_move_mean);

    % HFA instruct vs move
    subplot(params.fig.zcorr.subplot_x,params.fig.zcorr.subplot_y,1)
    x = nanmean(thisHFA_instruct(:,tStart_early:tEnd_late),2);
    y = nanmean(thisHFA_move(:,tStart_early:tEnd_late),2);
    scatter(x,y,50,'filled',...
        'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.3 .3 .3])
    [r,p] = corr (x,y,'rows','complete');
    [b] = regress(y,[ones(size(x,1),1) [x]]);
    [top_int, bot_int,X] = regression_line_ci(0.05,b,x,y);
    swagAxes(gca,params.fig.zcorr.fSize,'z HFA instruct','z HFA move',...
        ['corr r = ' num2str(r) '; p = ' num2str(p)]);
    
    % HFA move  vs. beta move 
    subplot(params.fig.zcorr.subplot_x,params.fig.zcorr.subplot_y,2)
    x = nanmean(thisHFA_move(:,tStart_early:tEnd_late),2);
    y = nanmean(thisBeta_move(:,tStart_early:tEnd_late),2);
    scatter(x,y,50,'filled',...
        'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.3 .3 .3])
    [r,p] = corr (x,y,'rows','complete');
    [b] = regress(y,[ones(size(x,1),1) [x]]);
    [top_int, bot_int,X] = regression_line_ci(0.05,b,x,y);
    swagAxes(gca,params.fig.zcorr.fSize,'z HFA move','z beta move',...
        ['corr r = ' num2str(r) '; p = ' num2str(p)]);

    
%     % HFA move  vs. beta move 
%     subplot(params.fig.zcorr.subplot_x,params.fig.zcorr.subplot_y,4)
%     x = nanmean(thisHFA_move(:,tStart_early:tEnd_late),2);
%     y = nanmean(thisBeta_move(:,tStart_early:tEnd_late),2);
%     scatter(x,y,50,'filled',...
%         'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.3 .3 .3])
%     [r,p] = corr (x,y,'rows','complete');
%     [b] = regress(y,[ones(size(x,1),1) [x]]);
%     [top_int, bot_int,X] = regression_line_ci(0.05,b,x,y);
%     swagAxes(gca,params.fig.zcorr.fSize,'z HFA move','z beta move',...
%         ['corr r = ' num2str(r) '; p = ' num2str(p)]);
%     
    
%     % Beta instruct vs. move
%     subplot(params.fig.zcorr.subplot_x,params.fig.zcorr.subplot_y,2)
%     x = nanmean(thisBeta_instruct(:,tStart_late:tEnd_late),2);
%     y = nanmean(thisBeta_move(:,tStart_late:tEnd_late),2);
%     scatter(x,y,50,'filled',...
%         'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.3 .3 .3])
%     [b] = regress(y,[ones(size(x,1),1) [x]]);
%     [top_int, bot_int,X] = regression_line_ci(0.05,b,x,y);
%     [r,p] = corr (x,y,'rows','complete');
%     swagAxes(gca,params.fig.zcorr.fSize,'z beta instruct late','z beta move late',...
%         ['corr r = ' num2str(r) '; p = ' num2str(p)]);
%     
%     
%     % Theta instruct early vs. move late
%     subplot(params.fig.zcorr.subplot_x,params.fig.zcorr.subplot_y,3)
%     x = nanmean(thisTheta_instruct(:,tStart_early:tEnd_early),2);
%     y = nanmean(thisTheta_move(:,tStart_late:tEnd_late),2);
%     scatter(x,y,50,'filled',...
%         'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.3 .3 .3])
%     [r,p] = corr (x,y,'rows','complete');
%     [b] = regress(y,[ones(size(x,1),1) [x]]);
%     [top_int, bot_int,X] = regression_line_ci(0.05,b,x,y);
%     swagAxes(gca,params.fig.zcorr.fSize,'z theta instruct early','z theta move late',...
%         ['corr r = ' num2str(r) '; p = ' num2str(p)]);
%     
%     % HFA move late vs. beta move late
%     subplot(params.fig.zcorr.subplot_x,params.fig.zcorr.subplot_y,4)
%     x = nanmean(thisHFA_move(:,tStart_late:tEnd_late),2);
%     y = nanmean(thisBeta_move(:,tStart_late:tEnd_late),2);
%     scatter(x,y,50,'filled',...
%         'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.3 .3 .3])
%     [r,p] = corr (x,y,'rows','complete');
%     [b] = regress(y,[ones(size(x,1),1) [x]]);
%     [top_int, bot_int,X] = regression_line_ci(0.05,b,x,y);
%     swagAxes(gca,params.fig.zcorr.fSize,'z HFA move late','z beta move late',...
%         ['corr r = ' num2str(r) '; p = ' num2str(p)]);
%     
%     % Theta move late vs. beta move late
%     subplot(params.fig.zcorr.subplot_x,params.fig.zcorr.subplot_y,5)
%     x = nanmean(thisTheta_move(:,tStart_late:tEnd_late),2);
%     y = nanmean(thisBeta_move(:,tStart_late:tEnd_late),2);
%     scatter(x,y,50,'filled',...
%         'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.3 .3 .3])
%     [r,p] = corr (x,y,'rows','complete');
%     [b] = regress(y,[ones(size(x,1),1) [x]]);
%     [top_int, bot_int,X] = regression_line_ci(0.05,b,x,y);
%     swagAxes(gca,params.fig.zcorr.fSize,'z theta move late','z beta move late',...
%         ['corr r = ' num2str(r) '; p = ' num2str(p)]);
%             
%        
%      % Theta instruct early vs. beta instruct early
%     subplot(params.fig.zcorr.subplot_x,params.fig.zcorr.subplot_y,6)
%     x = nanmean(thisTheta_instruct(:,tStart_late:tEnd_late),2);
%     y = nanmean(thisBeta_instruct(:,tStart_late:tEnd_late),2);
%     scatter(x,y,50,'filled',...
%         'markerfacecolor',[.7 .7 .7],'markeredgecolor',[.3 .3 .3])
%     [r,p] = corr (x,y,'rows','complete');
%     [b] = regress(y,[ones(size(x,1),1) [x]]);
%     [top_int, bot_int,X] = regression_line_ci(0.05,b,x,y);
%     swagAxes(gca,params.fig.zcorr.fSize,'z theta instruct late','z beta instruct late',...
%         ['corr r = ' num2str(r) '; p = ' num2str(p)]);
%         

end
function [h,statsStruct] = zPowPlot_local(roiStruct,params,statsStruct)
tBins = roiStruct.all.tStruct_iw_grp.tBins;
fBins = roiStruct.all.tStruct_iw_grp.fBins;
fSize = params.fig.tf.fSize;

% create figure for each freq band
h(1) = swagFig([-0.0469    0.4975    1.1125    0.3825]); hold all; 
h(1).Name = 'zPowHFA';
h(2) = swagFig([-0.0469    0.4975    1.1125    0.3825]); hold all; 
h(2).Name = 'zPowBeta';
h(3) = swagFig([-0.0469    0.4975    1.1125    0.3825]); hold all; 
h(3).Name = 'zPowTheta';



% loop through regions
for i = 1:length(params.filt_roi)
    % get vars for this ROI
    roi = params.filt_roiLbl{i};

    % skip this if no activations
    if ~isfield(roiStruct,params.filt_roiLbl{i})
        continue
    end
       
    % HFA 0 - 500 MW - IW
    thisHFA_wait = cat(1,roiStruct.(roi).tStruct_iw_grp.zHFA_wait_mean);
    thisHFA_instruct = cat(1,roiStruct.(roi).tStruct_iw_grp.zHFA_instruct_mean);
    thisHFA_move = cat(1,roiStruct.(roi).tStruct_iw_grp.zHFA_move_mean);

    thisBeta_wait = cat(1,roiStruct.(roi).tStruct_iw_grp.zBeta_wait_mean);
    thisBeta_instruct = cat(1,roiStruct.(roi).tStruct_iw_grp.zBeta_instruct_mean);
    thisBeta_move = cat(1,roiStruct.(roi).tStruct_iw_grp.zBeta_move_mean);

    thisTheta_wait = cat(1,roiStruct.(roi).tStruct_iw_grp.zTheta_wait_mean);
    thisTheta_instruct = cat(1,roiStruct.(roi).tStruct_iw_grp.zTheta_instruct_mean);
    thisTheta_move = cat(1,roiStruct.(roi).tStruct_iw_grp.zTheta_move_mean);
    

    
    % HFA
    thisH = h(1);
    thisPowLbl = 'HFA';
    zPow_wait = thisHFA_wait;
    zPow_instruct = thisHFA_instruct;
    zPow_move = thisHFA_move;
    [statsStruct,b(i)] = doZPowComparison_local(roiStruct.(roi),thisH, thisPowLbl,...
        zPow_wait,zPow_instruct,zPow_move,statsStruct,params,tBins,params.filt_roi_kol{i});

    % Beta
    thisH = h(2);
    thisPowLbl = 'beta';
    zPow_wait = thisBeta_wait;
    zPow_instruct = thisBeta_instruct;
    zPow_move = thisBeta_move;
    [statsStruct,bb(i)] = doZPowComparison_local(roiStruct.(roi),thisH, thisPowLbl,...
        zPow_wait,zPow_instruct,zPow_move,statsStruct,params,tBins,params.filt_roi_kol{i});

    % Theta
    thisH = h(3);
    thisPowLbl = 'theta';
    zPow_wait = thisTheta_wait;
    zPow_instruct = thisTheta_instruct;
    zPow_move = thisTheta_move;
    [statsStruct,bbb(i)] = doZPowComparison_local(roiStruct.(roi),thisH, thisPowLbl,...
        zPow_wait,zPow_instruct,zPow_move,statsStruct,params,tBins,params.filt_roi_kol{i});
       
end

%place legend
for i = 1:3
    figure(h(i)); 
    subplot(params.fig.zpow.subplot_x,params.fig.zpow.subplot_y,3); hold all
    l = legend(b,params.filt_roi,...
        'location','none','fontsize',params.fig.zpow.fSize,'position',[.95 .8 0 0]);

    
    l = legend(bb,params.filt_roi,...
        'location','none','fontsize',params.fig.zpow.fSize,'position',[.95 .8 0 0]);

    l = legend(bbb,params.filt_roi,...
        'location','none','fontsize',params.fig.zpow.fSize,'position',[.95 .8 0 0]);
end


function [statsStruct,thisB] = doZPowComparison_local(thisSpecStruct,thisH, thisPowLbl,...
        zPow_wait,zPow_instruct,zPow_move,statsStruct,params,tBins,kol)
% collapse within subj
if params.fig.zpow.collapseWinSubj ==1
    [zPow_wait] = collapseWinSubj_local(thisSpecStruct,zPow_wait,2);
    [zPow_instruct] = collapseWinSubj_local(thisSpecStruct,zPow_instruct,2);
    [zPow_move] = collapseWinSubj_local(thisSpecStruct,zPow_move,2);
end

% HFA
figure(thisH)

% WAIT 
subplot(params.fig.zpow.subplot_x,params.fig.zpow.subplot_y,1); hold all
b(1) = plotSEM(tBins,nanmean(zPow_wait,1),SEM(zPow_wait),kol);
swagAxes(gca,params.fig.zpow.fSize,'Time from WAIT cue (ms)',thisPowLbl)
set(gca,'ylim',params.fig.zpow.ylim)
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [.4 .4 .4],'linestyle','--')
plot(get(gca,'xlim'),[0 0],'linewidth',3,'color', [.4 .4 .4],'linestyle','--')

% INSTRUCT
subplot(params.fig.zpow.subplot_x,params.fig.zpow.subplot_y,2); hold all
b(2) = plotSEM(tBins,nanmean(zPow_instruct,1),SEM(zPow_instruct),kol);
swagAxes(gca,params.fig.zpow.fSize,'Time from INSTRUCT cue (ms)',thisPowLbl)
set(gca,'ylim',params.fig.zpow.ylim)
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [.4 .4 .4],'linestyle','--')
plot(get(gca,'xlim'),[0 0],'linewidth',3,'color', [.4 .4 .4],'linestyle','--')

% MOVE
subplot(params.fig.zpow.subplot_x,params.fig.zpow.subplot_y,3); hold all
b(3) = plotSEM(tBins,nanmean(zPow_move,1),SEM(zPow_move),kol);
swagAxes(gca,params.fig.zpow.fSize,'Time from MOVE cue (ms)',thisPowLbl)
set(gca,'ylim',params.fig.zpow.ylim)
% plot line
plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [.4 .4 .4],'linestyle','--')
plot(get(gca,'xlim'),[0 0],'linewidth',3,'color', [.4 .4 .4],'linestyle','--')


 for bb = 1:length(b)
    set(b(bb),'facealpha',.75)
 end
 
 thisB = b(end);

function [h,statsStruct] = tfPlot_local(specStruct_ROI,params)

tBins = specStruct_ROI.all.tStruct_iw_grp.tBins;
fBins = specStruct_ROI.all.tStruct_iw_grp.fBins;
fSize = params.fig.tf.fSize;
clims = params.fig.tf.clim; 

% loop through regions
for i = 1:length(params.filt_roi)
    % create figure for region
    h(i) = swagFig([-0.0469    0.4975    1.1125    0.3825]); hold all; 
    h(i).Name = ['TF-' params.filt_roiLbl{i}];
    % get vars for this ROI
    roi = params.filt_roiLbl{i};

    % skip this if no activations
    if ~isfield(specStruct_ROI,params.filt_roiLbl{i})
        continue
    end
   
    % get matrix of power differences for instruct - wait
    tMat_iw = cat(3,specStruct_ROI.(roi).tStruct_iw_grp.tMat);
    tMat_mw = cat(3,specStruct_ROI.(roi).tStruct_mw_grp.tMat);
    tMat_mi = cat(3,specStruct_ROI.(roi).tStruct_mi_grp.tMat);

%     tMat_iw = cat(3,specStruct_ROI.(roi).tStruct_iw_grp.zPowDiffMat);
%     tMat_mw = cat(3,specStruct_ROI.(roi).tStruct_mw_grp.zPowDiffMat);
%     tMat_mi = cat(3,specStruct_ROI.(roi).tStruct_mi_grp.zPowDiffMat);
% 
    % collapse w/in subj
    if params.fig.tf.collapseWinSubj ==1
        [tMat_iw] = collapseWinSubj_local(specStruct_ROI.(roi),tMat_iw,3);
        [tMat_mw] = collapseWinSubj_local(specStruct_ROI.(roi),tMat_mw,3);
        [tMat_mi] = collapseWinSubj_local(specStruct_ROI.(roi),tMat_mi,3);
    end
    
    % do cluster stats on group data
    [clusStats_iw] = doComparison_local(tMat_iw,zeros(size(tMat_iw)),params,fBins,tBins);
    [clusStats_mw] = doComparison_local(tMat_mw,zeros(size(tMat_mw)),params,fBins,tBins);
    [clusStats_mi] = doComparison_local(tMat_mi,zeros(size(tMat_mi)),params,fBins,tBins);
    
    % update statsStruct
    statsStruct.clusStats.(roi).IW = clusStats_iw;
    statsStruct.clusStats.(roi).MW = clusStats_mw;
    statsStruct.clusStats.(roi).MI = clusStats_mi;
    
    % Time freq plot -- IW
    compLbl = 'Instruct--Wait';
    figure(h(i))
    subplot(params.fig.tf.subplot_x,params.fig.tf.subplot_y,1); hold all
    nElecs = length(specStruct_ROI.(roi).tStruct_iw_grp);
    nSubj = length(unique({specStruct_ROI.(roi).tStruct_iw_grp.subj}));
    plotThisTF_local(tBins,fBins,clusStats_iw.tMat_sig,clims,compLbl,nElecs,nSubj,fSize,roi)
    
    % Time freq plot -- MW
    compLbl = 'Move--Wait';
    figure(h(i))
    subplot(params.fig.tf.subplot_x,params.fig.tf.subplot_y,2); hold all
    nElecs = length(specStruct_ROI.(roi).tStruct_mw_grp);
    nSubj = length(unique({specStruct_ROI.(roi).tStruct_mw_grp.subj}));
    plotThisTF_local(tBins,fBins,clusStats_mw.tMat_sig,clims,compLbl,nElecs,nSubj,fSize,roi)
    
     % Time freq plot -- MI
    compLbl = 'Move--Instruct';
    figure(h(i))
    subplot(params.fig.tf.subplot_x,params.fig.tf.subplot_y,3); hold all
    nElecs = length(specStruct_ROI.(roi).tStruct_mw_grp);
    nSubj = length(unique({specStruct_ROI.(roi).tStruct_mi_grp.subj}));
    plotThisTF_local(tBins,fBins,clusStats_mi.tMat_sig,clims,compLbl,nElecs,nSubj,fSize,roi)
end


function [] = plotThisTF_local(tBins,fBins,tMat,clims,compLbl,nElecs,nSubj,fSize,roi)
    imagescnan(tBins,1:length(fBins),tMat);                 
    set(gca,'ydir','normal','clim',clims)
    setFreqYTicks(gca,fBins)
    swagAxes(gca,fSize,'Time from cue (ms)','Frequency (Hz)',...
    [roi ' z Power ' compLbl '(# elec = ' num2str(nElecs) ', # subj = ' num2str(nSubj) ')'])     
    % plot line
    plot([0 0],get(gca,'ylim'),'linewidth',3,'color', [1 1 1])
    c = colorbar('eastoutside');colormap('jet')
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


function [numSigEffects,matSigEffects] = getCoactivationData_local(thisClusStruct,params)
freqBandLbls = fieldnames(params.freqBandLims);
% get list of unique electrodes
uElbl_unique = unique({thisClusStruct.uElbl});

% make a matrix of
matSigEffects = zeros(length(uElbl_unique),length(freqBandLbls));

% number of activations
numSigEffects = nan(1,length(uElbl_unique));


for i = 1:length(uElbl_unique)
    uIdx = strcmp(uElbl_unique{i},{thisClusStruct.uElbl});
    thisClusStruct_e = thisClusStruct(uIdx);
    numSigEffects(i) = length(thisClusStruct_e);
    for e = 1:length(thisClusStruct_e)
        % populate matrix 
        fIdx = (find(strcmp(thisClusStruct_e(e).freqBandLbl,freqBandLbls)));
        matSigEffects(i,fIdx) = thisClusStruct_e(e).clusPolarity;
           
    end
end

function [h_coactivation,h_freqIncDec,h_sel,h_timing_HFA] = freqBandPlots_local(clusStruct_ROI,params,nsubj)

if params.fig.counts.subplot_x<params.fig.counts.subplot_y
    h_coactivation =  swagFig([-0.0469    0.4975    1.1125    0.3825]); hold all; 
    h_freqIncDec =  swagFig([-0.0469    0.4975    1.1125    0.3825]); hold all; 
    h_sel = swagFig([-0.0469    0.4975    1.1125    0.3825]); hold all; 
    h_zPowByBand = swagFig([-0.0469    0.4975    1.1125    0.3825]); hold all; 
else
    h_coactivation =  swagFig([0    0.0638    0.7164    0.8163]); hold all; 
    h_freqIncDec =  swagFig([0    0.0638    0.7164    0.8163]); hold all; 
    h_sel = swagFig([0    0.0638    0.7164    0.8163]); hold all;     
    h_zPowByBand = swagFig([0    0.0638    0.7164    0.8163]); hold all;
end

h_timing_HFA = swagFig([0.6617    0.4300    0.3250    0.4500]); hold all; 

% make a freq span and mean freq matrix for freq bands of interest
fBandLbls = fieldnames(params.freqBandLims);
tBins = clusStruct_ROI.all.tBins;
fBins = clusStruct_ROI.all.fBins;
compLbl = [clusStruct_ROI.all(1).retLbl1 '-' clusStruct_ROI.all(1).retLbl2];

for i = 1:length(params.filt_roi)    
    % get vars for this ROI
    roi = params.filt_roiLbl{i};
    
    % skip this if no activations
    if ~isfield(clusStruct_ROI,roi)
        continue
    end
   
    % get co-activation data (per electrode)
    [numSigEffects.(roi),matSigEffects.(roi)] = ...
        getCoactivationData_local(clusStruct_ROI.(roi),params);
    
    % variables used for counts analysis
    isPos = [clusStruct_ROI.(roi).clusPolarity]'==1&[clusStruct_ROI.(roi).p]'<=params.pThresh; % increase 
    isNeg = [clusStruct_ROI.(roi).clusPolarity]'==-1&[clusStruct_ROI.(roi).p]'<=params.pThresh;%vs. decrease
    isSelective = [clusStruct_ROI.(roi).anova_dirSelectivity_p]'<=params.pThresh_moveSelective;    
    zPowDiffVecMean = cat(1,clusStruct_ROI.(roi).zPowDiffVecMean);
    zPowDiffVecMean(isinf(zPowDiffVecMean)) = nan;
    zPowDiff = cat(1,clusStruct_ROI.(roi).zPowDiffClusMean);
    zPowDiff(isinf(zPowDiff)) = nan;
    % do counts
    for f = 1:length(fBandLbls)
        % assigns each signficant effect to a frequency band 
        retIdx = strcmp({clusStruct_ROI.(roi).freqBandLbl},fBandLbls{f})';
       
        % freq inc/dec
        freqBandIncreases_prctSig.(roi).(fBandLbls{f}) = (sum(retIdx&isPos)/length(retIdx))*100;       
        freqBandDecreases_prctSig.(roi).(fBandLbls{f}) = (sum(retIdx&isNeg)/length(retIdx))*100;
        
        % selectivity
        freqBandSelective_prctSig.(roi).(fBandLbls{f}) = (sum(retIdx&isSelective)/sum(retIdx))*100;
   
        % mean timing
        thisPow = zPowDiff(retIdx);
        thisPowVec = zPowDiffVecMean(retIdx,:);       
        if params.fig.timing.collapseWinSubj ==1
          [thisPow] = collapseWinSubj_local(clusStruct_ROI.(roi)(retIdx),thisPow,2);
          [thisPowVec] = collapseWinSubj_local(clusStruct_ROI.(roi)(retIdx),thisPowVec,2);
        end
        freqBandPow.(roi).(fBandLbls{f}) = thisPow;
        freqBandTimeCourse.(roi).(fBandLbls{f}) = thisPowVec;       
    end
       
    % set up vars for bar plot
    thisCount_sel = struct2cell(freqBandSelective_prctSig.(roi));
    thisCount_sel= [thisCount_sel{:}];
    
    thisCount_pos = struct2cell(freqBandIncreases_prctSig.(roi));
    thisCount_pos= [thisCount_pos{:}];
    
    thisCount_neg = struct2cell(freqBandDecreases_prctSig.(roi));
    thisCount_neg= [thisCount_neg{:}];   
    
    % plot co-activation data
    figure(h_coactivation)
    subplot(params.fig.counts.subplot_x,params.fig.counts.subplot_y,i); hold all
    counts = (histcounts(numSigEffects.(roi))/length(numSigEffects.(roi)))*100;
    b = bar(counts);
    set(b(1),'facecolor',[.7 .7 .7],'edgecolor',[.2 .2 .2],'linewidth',2)
    set(gca,'xtick',[1:length(counts)]);
    swagAxes(gca,params.fig.counts.fSize,'# sig effects per electrode','% all significant electrodes within region',...
    [roi '(# sig electrodes = ' num2str(length(numSigEffects.(roi))) ', # subj = ' num2str(nsubj.(roi)) ')'])     
    set(gca,'xlim',[.5 length(counts)+.5],'ylim',[0 70])
  
    %subplot frequency band counts
    figure(h_freqIncDec)
    subplot(params.fig.counts.subplot_x,params.fig.counts.subplot_y,i); hold all
    b = bar([thisCount_pos;thisCount_neg]','stacked');
    set(b(1),'facecolor',[1 .7 .7])
    set(b(2),'facecolor',[.7 .7 1])
    set(gca,'xtick',[1:length(fBandLbls)],'xticklabel',fBandLbls);
    swagAxes(gca,params.fig.counts.fSize,'Frequency band','% all significant effects within region',...
    [roi '(# sig effects = ' num2str(length(retIdx)) ', # subj = ' num2str(nsubj.(roi)) ')'])     
    set(gca,'xlim',[.5 length(fBandLbls)+.5],'ylim',[0 50])
    legend(b, {'Power increases','Power decreases'})
    
    figure(h_sel)
    subplot(params.fig.counts.subplot_x,params.fig.counts.subplot_y,i); hold all
    b = bar(thisCount_sel);
    set(b(1),'edgecolor','k','facecolor',[.7 .7 .7],'linewidth',1.5)
    set(gca,'xtick',[1:length(fBandLbls)],'xticklabel',fBandLbls);
    swagAxes(gca,params.fig.counts.fSize,'Frequency band','% movement selective within frequency and region',...
    [roi '(# sig effects = ' num2str(length(retIdx)) ', # subj = ' num2str(nsubj.(roi)) ')'])     
    set(gca,'xlim',[.5 length(fBandLbls)+.5],'ylim',[0 70])

    % plot timing
    figure(h_timing_HFA);hold all
    thisFbandLbl = 'HFA'; 
    powVec = abs(freqBandTimeCourse.(roi).(thisFbandLbl)); 
    plot(tBins',nanmean(powVec,1),'--k','linewidth',2)
    b_timing_HFA(i) = plotSEM(tBins',nanmean(powVec,1),SEM(powVec),params.fig.timing.cmap(i,:));
    set(b_timing_HFA(i),'facealpha',0.6)
    swagAxes(gca,params.fig.timing.fSize,'Time from cue',['z ' thisFbandLbl ' abs(' compLbl ')'])
    plot([0 0],get(gca,'ylim'),'linestyle','--','linewidth',3,'color', [.5 .5 .5])

    % plot mean z-pow within each band
    figure(h_zPowByBand)
    subplot(params.fig.counts.subplot_x,params.fig.counts.subplot_y,i); hold all
    thisPowDiff = struct2cell(freqBandPow.(roi));
    thisPowDiff_mean = cellfun(@nanmean,thisPowDiff);
    thisPowDiff_sem = cellfun(@SEM,thisPowDiff);
    b = bar(thisPowDiff_mean);
    e = errorbar(thisPowDiff_mean,thisPowDiff_sem,'.k','linewidth',1.5);
    set(b(1),'edgecolor','k','facecolor',[.7 .7 .7],'linewidth',1.5)
    set(gca,'xtick',[1:length(fBandLbls)],'xticklabel',fBandLbls);
    swagAxes(gca,params.fig.counts.fSize,'Frequency band','Mean z power difference',...
    [roi '(# sig effects = ' num2str(length(retIdx)) ', # subj = ' num2str(nsubj.(roi)) ')'])     
    set(gca,'xlim',[.5 length(fBandLbls)+.5])

end

% % set legend for timing figure
figure(h_timing_HFA)    
legend(b_timing_HFA,params.filt_roiLbl,'location','northwest')
 
function [thisVar_subj] = collapseWinSubj_local(thisSpecStruct,thisVar,dim)
% dim is the non-subject dimension 
subjList_unique = unique({thisSpecStruct.tStruct_iw_grp.subj});
subjList_all = {thisSpecStruct.tStruct_iw_grp.subj};
switch dim
    case 1
        thisVar_subj = nan(1,length(subjList_unique));
    case 2
        thisVar_subj = nan(length(subjList_unique),size(thisVar,2));
    case 3
        thisVar_subj = nan(size(thisVar,1),size(thisVar,2),length(subjList_unique));


end
for s = 1:length(subjList_unique)
    sIdx = strcmp(subjList_unique{s},subjList_all);
    switch dim
        case 1
            thisVar_subj(s) = nanmean(thisVar(sIdx));
        case 2
            thisVar_subj(s,:) = nanmean(thisVar(sIdx,:),1);
        case 3
            thisVar_subj(:,:,s) = nanmean(thisVar(:,:,sIdx),3);
    end
end


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


function [] = setFreqYTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'ytick',yt,...
'yticklabel',round(fBins(yt)))


function [] = setFreqXTicks(ax,fBins)
yt = [[1:5:length(fBins)] length(fBins)];
if any(diff(yt)==0)
    yt = [1:5:length(fBins)];
end
set(ax,'xtick',yt,...
'xticklabel',round(fBins(yt)))
h=fill([x,flipdim(x,2)],[y+SEM,flipdim(y-SEM,2)]);
%% Cluster based statistics subfunctions

function [truStats] = doComparison_local(powVec1,powVec2,params,fBins,tBins);

% get clus alpha values
clusConfig = params.clusConfig;

[truStats.tMat,truStats.pMat,truStats.df,...
    truStats.posClusStats,truStats.negClusStats,...
    truStats.posClusInd,truStats.negClusInd] = ...
    clusStats_wrapper_local(powVec1,powVec2,params);
   
% shufled power
shufPow = cat(3,powVec1,powVec2);

for i = 1:clusConfig.iters
    %shuf fb events
    shufInd = logical(randi([0 1],[1 size(shufPow,3)]));
    
    [~,~,~,shuf_posClusStats,shuf_negClusStats] = ...
    clusStats_wrapper_local(shufPow(:,:,shufInd),shufPow(:,:,~shufInd),params);
    
    shufStats(i).posClusStats = shuf_posClusStats(1);
    shufStats(i).negClusStats = shuf_negClusStats(1);
    
    if mod(i,100) == 0
      disp(i)
    end
end


% create tmat sig (for plotting)
truStats.tMat_sig = nan(size(truStats.tMat));


% calc p-values and adjust indices (and update tMat_sig)
%initialize
truStats.posClusStatsP = [];truStats.negClusStatsP = [];
truStats.posClusFreqBandLbl = cell(size(truStats.posClusStatsP));
truStats.negClusFreqBandLbl = cell(size(truStats.negClusStatsP));
truStats.posClusTimeBandLbl = cell(size(truStats.posClusStatsP));
truStats.negClusTimeBandLbl = cell(size(truStats.negClusStatsP));

for t = 1:length(truStats.posClusStats)
    [truStats.posClusStatsP(t)] = calcP(truStats.posClusStats(t),[shufStats.posClusStats]);
    [truStats.posClusInd{t}] = truStats.posClusInd{t};
    if truStats.posClusStatsP(t) <= (clusConfig.alfa(1))*2
        truStats.tMat_sig(truStats.posClusInd{t}) = truStats.tMat(truStats.posClusInd{t});
        [truStats] = assignFreqTimBandLbls_local(truStats,'pos',t,fBins,tBins,params);
    end
end
for t = 1:length(truStats.negClusStats)
    [truStats.negClusStatsP(t)] = calcP(truStats.negClusStats(t),[shufStats.negClusStats]);
    [truStats.negClusInd{t}] =  truStats.negClusInd{t};
    if truStats.negClusStatsP(t) <= (clusConfig.alfa(1))*2
        truStats.tMat_sig(truStats.negClusInd{t}) = truStats.tMat(truStats.negClusInd{t});
        [truStats] = assignFreqTimBandLbls_local(truStats,'neg',t,fBins,tBins,params);

    end
end

%update truStats with shufStats so we can plot them later to show false pos
%rates
truStats.shufStats_pos = [shufStats.posClusStats];
truStats.shufStats_neg = [shufStats.negClusStats];


function [truStats] = assignFreqTimBandLbls_local(truStats,posNegStr,c,fBins,tBins,params);

switch posNegStr
    case 'pos'
       thisClusInd = truStats.posClusInd{c}; 
    case 'neg'
        thisClusInd = truStats.negClusInd{c}; 
end


[truStats.freqIdx,truStats.timeIdx] = ...
    ind2sub(size(truStats.tMat),thisClusInd);

truStats.freqs = fBins(truStats.freqIdx)';
truStats.times = tBins(truStats.timeIdx);
truStats.freqMax = max(truStats.freqs);
truStats.freqMin = min(truStats.freqs);
truStats.freqSpan = max(truStats.freqs)-min(truStats.freqs);
truStats.timeMax = max(truStats.times);
truStats.timeMin = min(truStats.times);



% % assign a frequncy band and time band label and update tStruct as well
fbandLbls = fieldnames(params.freqBandLims);
fRange = [truStats.freqMin truStats.freqMax];
freqBandLims_mat = struct2cell(params.freqBandLims);
freqBandLims_mat = cat(1,freqBandLims_mat{:});

tBandLbls = fieldnames(params.timeBandlims);
tRange = [truStats.timeMin truStats.timeMax];
timeBandLims_mat = struct2cell(params.timeBandlims);
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

switch  posNegStr
    case 'pos'
       truStats.posClusFreqBandLbl{c} = fbandLbls{fBandIdx};
       truStats.posClusTimeBandLbl{c} = tBandLbls{tBandIdx};
       
    case 'neg'
       truStats.negClusFreqBandLbl{c} = fbandLbls{fBandIdx};
       truStats.negClusTimeBandLbl{c} = tBandLbls{tBandIdx};
end

function [tVec,pVec,df,...
    posClusStats,negClusStats,posClusInd,negClusInd] = clusStats_wrapper_local(powVec1,powVec2,params)
% get clus alpha values
clusConfig = params.clusConfig;

% do one-tailed t-test
[~,pVec,~,stats] = ttest2(powVec1,powVec2,[],'right',[],3);
tVec = stats.tstat;
df = stats.df;




%calc clus stat
[posClusStats, negClusStats, posClusInd, negClusInd] = ...
    calcClusStat(tVec,pVec,clusConfig.alfa);



function[posClusStats, negClusStats, posClusInd, negClusInd] = calcClusStat(tMat,pMat,alfa)
%This function identifies contiguous tf clusters. 
%Inputs:
%tMat: tf matrix of t values
%pMat = tf matrix of p values (one-tailed)
%alfaVec = [0.025 0.975]; one-tailed p-thresholds; p<0.025  will identify
           %positive tf activations,  p>0.975 will reveal negative tf activations
%Outputs:
%posClusStats .... vector of positive activations ordered by magniture
%negClusStats ..... vector of negative activations ordered by magnitude
%posClusInd ...... cell array of sig. indices for positive activations
%negClusInd ..... cell array of sig. indices for negative activaions 

%Positive activation
%identify sig. windows (e.g., p < 0.05)
inds = find(pMat<alfa(1));
[posClusStats,posClusInd] = doClusCalc(inds,tMat,pMat);

inds = find(pMat>alfa(2));
[negClusStats,negClusInd] = doClusCalc(inds,tMat,pMat);

function[clusStats,clusInds] = doClusCalc(inds,tMat,pMat)
if length(inds)<2 % if 
    clusStats = [nan];
    clusInds = {[]};
    return
end
%identify contiguous windows
[subs(:,1), subs(:,2)] = ind2sub(size(pMat),inds); %this step will be helpful for spectral analyses
z = linkage(subs);
clusObs=cluster(z,'cutoff',sqrt(2),'criterion','distance');

%calc time of onset and cluster stats for each ''activation''
clusStats = nan(1,length(unique(clusObs)));
clusInds = cell(1,length(unique(clusObs)));
for i = 1:length(unique(clusObs))
    clusStats(i) = sum(tMat(inds(clusObs==i)));
    clusInds{i} = inds(clusObs==i);
end

%sort by magnitude
[~,magInd] = sort(abs(clusStats),'descend');
clusStats = clusStats(magInd);
clusInds = clusInds(magInd);


function [pVal] = calcP(true_t,shuf_t)
%For t-stats, it is a one-tailed p value; (high t gives low p, low t gives
%high p) 
%For Cluster t-stats, it is a two-tailed p-value for pos. clusters and 
%For f-stats it is a two-tailed pvalue; low ps are sig.
numDim = sum(size(true_t)>1);
if numDim == 0 %a single distribution of t-values (e.g clus stats)
    shuf_t = squeeze(shuf_t);
    if true_t >= 0 % for pos cluster stats
        pVal = sum(true_t<=shuf_t)./length(shuf_t);
    elseif true_t < 0 % for neg cluster stats
        pVal = sum(true_t>=shuf_t)./length(shuf_t);
    elseif isnan(true_t)
        pVal = nan;
        return
    end
else
    if numDim == 1
       tMat = repmat(true_t,[1 size(shuf_t,2) size(shuf_t,3)]);
    elseif numDim == 2
       tMat = repmat(true_t,[1 1 size(shuf_t,3)]);
    end
    %get a single null value at each freq 
    shuf_t = repmat(nanmean(shuf_t,2),[1 size(tMat,2) 1]);
    pVal = sum(tMat<=shuf_t,3)./size(shuf_t,3);
end

function [pStruct] = loadPStruct_local(params)

cd_mkdir(params.statsDir)
if exist('pStruct.mat','file')
    load('pStruct.mat')  
else
    pStruct = struct;
end

% function[clusConfig] = clusConfig_local();
% clusConfig.alfa = [0.025 0.975];
% clusConfig.iters = 1000;
%clusConfig.saveDir = fullfile(dirs.analyses,'ecog','subjPow','clusStats');
