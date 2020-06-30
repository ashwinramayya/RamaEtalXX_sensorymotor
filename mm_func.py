#import packages
import numpy as np # numpy
import pandas as pd # pandas
from scipy.io import loadmat # to load matlab
from itertools import compress
from scipy import stats
import mne
import pycircstat as circ
from mne import time_frequency
from fooof import FOOOF
from pylab import *
import os





######### Load Data ###########

def loadEEG(clusId=3):
    """ This function loads .mat files that were exported via matlab scripts that contain EEG traces that are organized into movement, wait, and instruction trials. They are pre-selected based on significant task-related power changes and organized into clusters based on spectro-temporal dynamics of those power changes. Cluster 1 shows instrcution-related HFA changes, Cluster 2 shows predominant theta related increases, and Cluster 3 shows predominant movement-related HFA increases. These data also contain information about subject, session, samplerate, and contact information. Line noise bandpass filters have already been applied at 60, 120 and 180 Hz

    #Parameters #clusId ..... Optional parameter. Default set to 3( can  also input 1,2)

    """

    # import packages
    #import numpy as np # numpy
    #import pandas as pd # pandas
    #from scipy.io import loadmat # to load matlab

    # create a params dictionary
    params = {'clusId_list':[1,2,3],'clusId':clusId,'contactIdx':0,
          'rawDir':'/Users/ashwinramayya/Dropbox/neurosurgery_research/data/scratch/clusStruct/pac/'}
    print(params['rawDir'])

    # Select FileName for selected cluster (#fname = 'pacClusStruct-3.mat' vs. 1 vs 2)
    fname = 'pacClusStruct-'+ str(params['clusId']) +'.mat'
    fullfile = params['rawDir']+fname;
    print(fullfile)


    # load matlab struct files (squeeze me returns a single dim array, rather than 2D)
    clusMat = loadmat(fullfile,squeeze_me=True)
    #see fields of loaded mat file
    print(clusMat.keys())

    # This data structure has n entrys (each representing a contact)
    num_contacts = clusMat['eegStruct_thisClus'].shape[0]
    print('Number of contacts in Cluster '+str(params['clusId'])+' = '+str(num_contacts))

    # Use the dtype field to see structures of the data struct (e.g.)
    #print(clusMat['eegStruct_thisClus'][0].dtype)


    # convert EEG data to ndarrays

    # create an empty list container to collect EEG data for all contacts in this cluster
    EEG = []

    # create an empty list contained to collect EEG metadata for each contact
    metadata = []#pd.DataFrame(data=None, index=None, columns=None, dtype=None)


    #list of fields to skip when collecting metadata
    exclude_list = ['eegConfig', 'EEG_wait',
                    'EEG_inst', 'EEG_move', 'rawEEG_wait', 'rawEEG_inst',
                    'rawEEG_move','tBins','fBins','zHFA_wait_mean','zHFA_instruct_mean',
                   'zHFA_move_mean','zBeta_wait_mean','zBeta_instruct_mean','zBeta_move_mean',
                   'zTheta_wait_mean','zTheta_instruct_mean','zTheta_move_mean']


    # Loop through contacts
    for c in iter(np.arange(num_contacts)):

        # check that loop is working correctly
        params['contactIdx'] = c
        #print(c)

        # select data from one contact
        thisContactData = clusMat['eegStruct_thisClus'][c]

        #  collect EEG data from this contact in an array
        thisContactEEG = np.array([(thisContactData['EEG_wait']),
                                   (thisContactData['EEG_inst']),
                                   (thisContactData['EEG_move'])],
                                  dtype='float')

        # add data to EEG container
        EEG.append(thisContactEEG)

        # create empty dictionary for eeg meta data
        eeg_meta_dict=({'':''})

        for i in thisContactData.dtype.names:
            if np.shape(thisContactData[i]) == (0,):        # replace empty vales with none
                print('REPLACING Contact '+str(c)+' field '+str(i)+' with NONE')
                eeg_meta_dict.update({i:'None'})
            # this deals with SR that are sometimes arrays
            elif i == 'SR' and (isinstance(thisContactData[i], int)==0):
                eeg_meta_dict.update({i:thisContactData[i][0]})
            elif (i in exclude_list)==0:
                eeg_meta_dict.update({i:thisContactData[i]})

        for i in thisContactData['eegConfig'].dtype.names:
            #print(i)
            try:
                eeg_meta_dict.update({i:thisContactData['eegConfig'][i].astype(np.int)})
            except:
                eeg_meta_dict.update({i:thisContactData['eegConfig'][i]})

        # EEG metadata to dataframe
        metadata.append(pd.DataFrame(eeg_meta_dict,index=[c]))

    # convert list to dataframe
    metadata = pd.concat(metadata)

    #
    print('Done loading Cluster '+str(params['clusId']))
    return [EEG,metadata]

def applyMasks(EEG,metadata,subj_mask=None,eLbl_mask=None,roi_mask=None):
    """Returns a list of contacts based on the masks that are provided as inputs"""

    # import packages
    #import numpy as np # numpy
    #import pandas as pd # pandas
    #from itertools import compress

    # create a params dictionary that will hold the various masks that are passed in
    # Default values:
    #subj_mask = None
    #eLbl_mask = None
    #roi_mask = 'perirolandic'#None
    masks_dict = {'subj':subj_mask,'eLbl':eLbl_mask,'ROI':roi_mask}


    # identify an index that will include all elecs
    idx = np.ones(metadata.shape[0],dtype='bool')
    print(idx.sum())

    masks = {'subj':idx,'eLbl':idx,'ROI':idx}
    print(masks['subj'].sum())


    # create a mask for each of the fields above
    for i in masks_dict.keys():
        # if mask label is "None", fill mask dict with idx to include all contacts
        if masks_dict[i] != None:
            masks[i] = metadata[i] == masks_dict[i]

    #combine masks
    masks_combined = masks['subj'] & masks['eLbl'] & masks['ROI']

    # Mask the EEG list and the meta
    EEG = list(compress(EEG,masks_combined))
    metadata = metadata[masks_combined]

    print(len(EEG))
    return [EEG,metadata]

def getEEGClip(eeg_contact,metadata_contact,event_type,time_window_ms,buffer_ms = 0):
    """Input EEG and metadata from a single contact,event type ( wait vs. inst vs move) and time window (relative to cue onset) Returns a ntrial x time EEG array. Returns a list of contacts based on the masks that are provided as an input"""

    #example input values:
    #event_type = 'move'
    #time_window_ms = (-1000,1000)
    #buffer_ms = 0

    #filter eeg contact by event type and time window
    if event_type == 'wait':
        eeg_contact = eeg_contact[0]
    elif event_type == 'inst':
        eeg_contact = eeg_contact[1]
    elif event_type == 'move':
        eeg_contact = eeg_contact[2]
    else:
        print('ERROR: Incorrect event type info provided')
        #return None

    # get samplerate
    SR = metadata_contact['SR'].to_numpy().astype('float')

    # total time window of EEG_contact
    duration_samp= metadata_contact['durationMS'].values.astype('float')*(SR/1000)
    offset_samp= metadata_contact['offsetMS'].values.astype('float')*(SR/1000)
    orig_samp = np.arange(0,duration_samp)+offset_samp

    # Our time windows of interest in samples
    time_window_samp = time_window_ms*(SR/1000)

    # buffer in samples
    buffer_samp = np.ceil(buffer_ms*(SR/1000))

    # find start and end in samples (including buffer)
    samp_start = np.ceil(time_window_samp[0]-buffer_samp)
    samp_end = np.ceil(time_window_samp[1]+buffer_samp)


    #find clip start and end in samples as follows
    # find the absolute value of the difference
    #between the samp start were looking for and the original samples using
    # np.absolute, then finding the index of the miniminum value(0) using argmin

    clip_start_samp = np.absolute(orig_samp-samp_start).argmin()
    clip_end_samp = np.absolute(orig_samp-samp_end).argmin()

    # mask eeg_contact
    eeg_clip = eeg_contact[:,np.arange(clip_start_samp,clip_end_samp)]

    # remove trials with any nans (this is problematic for mne functions downstream)
    # find trials with nan
    idx_to_rem = np.any(np.isnan(eeg_clip),axis=1)

    # remove trials with nans
    eeg_clip = eeg_clip[idx_to_rem==False][:]

    # return
    return [eeg_clip,time_window_samp,SR]

def to_mne(eeg_clip,metadata_contact,time_window_ms,reject_thresh=0):
    """This function takes an eeg clip from the function above and converts it to an mne epochs array"""
    import mne

    # format data as (n_epochs, n_channels, n_times)
    n_epochs = np.shape(eeg_clip)[0]
    n_channels = 1
    n_times =  np.shape(eeg_clip)[1]
    data = np.reshape(eeg_clip,(n_epochs,n_channels,n_times))

    # create_info
    #get ch name, SR
    ch_name = [metadata_contact['uElbl'].to_numpy().astype('str')[0]]
    SR = metadata_contact['SR'].to_numpy().astype('float')
    info = mne.create_info(ch_name,SR,ch_types='ecog')

    # create mne epochs array (tmin is the offset in seconds)
    eeg_clip_epochs = mne.EpochsArray(data,info,events=None,tmin=time_window_ms[0]/1000,reject= dict(ecog = reject_thresh))

    # return
    return eeg_clip_epochs

# define wrapper function to return a eeg clip from the EEG_masked list using the params below
def getEEGClip_wrapper(eeg_contact,metadata_contact,event_type='move',time_window_ms=(-1000,200),rejection_thresh_sd = 20,buffer_ms = 0):
    """This is a wrapper function for get EEG clip and EEG clip to MNE """

    # inputs:
    #eeg_contact - filtered for a single contact (e.g., EEG_masked[c], where c is the index of the contact)
    #metadata_contact - filtered for a single contact (e.g., metadata_masked[c:c+1], where c is the index of the
    #e.g.,
    #eeg_contact = EEG_masked[c];
    #metadata_contact = metadata_masked[c:c+1];
    #event_type='move'
    #time_window_ms=(-1000,200)
    #reject_thresh_sd = 20; # This sets the number of standard deviations above the mean to reject



    # find rejection threshold value based on percentile of data across the entire contact
    reject_thresh = np.nanmean(eeg_contact) + (np.nanstd(eeg_contact)*rejection_thresh_sd)

    # get eeg clip
    [eeg_clip,time_window_samp,SR] = getEEGClip(eeg_contact,metadata_contact,event_type,time_window_ms,buffer_ms)

    # To MNE
    eeg_clip_epochs = to_mne(eeg_clip,metadata_contact,time_window_ms,reject_thresh)

    return eeg_clip_epochs



######### Electrode Class  ###########


# This is the workhorse object for analyses. It creates an instance of an electrode contact and can perform ERP, FOOOF, PSD, power, across trial phase locking (within and vs. other electrodes). It is flexible in event type (move vs. instruction) and time window. It stores outputs of above analyses and updates these results as new analyses are run.

# define class
class Electrode:
    """This is a class that collects EEG and metadata for a single contact and applies various functions
    for feature extraction and basic plotting """
    # Initializer / Instance Attributes
    def __init__(self,eeg,metadata):
        self.eeg = eeg
        self.metadata = metadata
        self.subj = metadata['subj'].iloc[0]
        self.eLbl = metadata['eLbl'].iloc[0]
        self.uElbl = metadata['uElbl'].iloc[0]
        self.roi = metadata['ROI'].iloc[0]

        #containers
        # params
        self.event_type = []
        self.time_window_ms =[]
        self.buffer_ms = []
        self.rejection_thresh_sd =[]
        self.frange = []

        # data (eeg, PSD)
        self.freqs = []
        self.eeg_clip = []
        self.psd_mean = []
        self.psd_log_mean = []
        self.psd_log_sem = []

        # fooof results (initialize parameter containers as nans in case there are no peaks)
        self.fm = [] # for foooof object
        self.aperiodic = np.empty([3])
        self.aperiodic[:] = nan
        self.theta =  self.aperiodic.copy()
        self.alpha =  self.aperiodic.copy()
        self.beta =  self.aperiodic.copy()
        self.gamma = self.aperiodic.copy()
        self.r_sq =  []
        self.fit_error = []


        # time-freq results
        self.freqs_wave = []
        self.time_window_samp = [] # this tracks time_window in samples assoc w power/phase matrix after buffer has been removed

        # Power
        self.z_powmat = []
        self.z_hfa = []
        self.z_theta = []
        self.z_alpha = []
        self.z_beta = []
        self.z_gamma = []

        #Phase
        self.all_phase = []
        self.all_phase_diff = []

        #PLV
        self.PLV_freq_interest = []
        self.PLV_time_win_interest = []
        self.PLV_p = []
        self.PLV_z = []
        self.PLV_rvl = []
        self.phs_avg = []


    # get eeg clip
    def getEEG(self,event_type=None,time_window_ms=None,buffer_ms=None,rejection_thresh_sd=None):

        if event_type == None:
            event_type = self.event_type
        if time_window_ms == None:
            time_window_ms = self.time_window_ms
        if buffer_ms == None:
            buffer_ms = self.buffer_ms
        if rejection_thresh_sd == None:
            rejection_thresh_sd = self.rejection_thresh_sd
        self.eeg_clip = getEEGClip_wrapper(self.eeg,self.metadata,event_type,time_window_ms,rejection_thresh_sd,buffer_ms)


        # update self
        self.event_type = event_type
        self.time_window_ms = time_window_ms
        self.buffer_ms = buffer_ms
        self.rejection_thresh_sd = rejection_thresh_sd

    # plot ERP
    def plotSingleTrial(self,trialnum=10,ax=None):
        """ inputs
        ax = axes to plot erp in """
        if ax==None:
            fig = figure()
            ax = subplot(111)
        eeg_data = self.eeg_clip.get_data()[trialnum].squeeze()
        xval_s = np.linspace(self.time_window_ms[0],self.time_window_ms[1],np.shape(eeg_data)[0])
        plot(xval_s,eeg_data)
        xlabel('Time from '+self.event_type+' (s)')
        ylabel('Voltage (uV)')
        axis('tight')
        vlines(0,ax.get_ylim()[0],ax.get_ylim()[1]);
        hlines(0,ax.get_xlim()[0],ax.get_xlim()[1]);

    # plot ERP
    def plotERP(self,ax=None):
        """ inputs
        ax = axes to plot erp in """
        if ax==None:
            fig = figure()
            ax = subplot(111)
        evoked = self.eeg_clip.average(picks=['ecog']);
        evoked.plot(axes=ax,picks=None,show=False,time_unit='s',hline=[0.0],units = dict(ecog='z-score'),
               scalings = dict(ecog=1));
        vlines(0,ax.get_ylim()[0],ax.get_ylim()[1]);
        ax.set_title('')

    # PSD
    def calcPSD(self,frange):
        psd,self.freqs = mne.time_frequency.psd_welch(self.eeg_clip,fmin=frange[0],fmax=frange[1])

        # calculate mean psd across trials for fooof (without logging, fooof applies log internally)
        self.psd_mean = psd.mean(0).squeeze()

        # calc log mean and log sem for plotting
        self.psd_log_mean = np.log(psd).mean(0).squeeze()
        self.psd_log_sem = stats.sem(np.log(psd),axis=0).squeeze()

        #
        self.frange = frange


    def plotPSD(self,ax=None):
        if ax==None:
            fig = figure()
            ax = subplot(111)
        plot(self.freqs,self.psd_log_mean,linewidth=3)
        fill_between(self.freqs,self.psd_log_mean+self.psd_log_sem,self.psd_log_mean-self.psd_log_sem,alpha=0.4)
        xlabel('Frequency (Hz)')
        ylabel('log power')

    # FOOOF
    def clearFOOOF(self):
        """This function clears prevously stored fooof results so that there is
        no mixing results when you perform a new fit """
        self.fm = [] # for foooof object
        self.aperiodic = np.empty([3])
        self.aperiodic[:] = nan
        self.theta =  self.aperiodic.copy()
        self.alpha =  self.aperiodic.copy()
        self.beta =  self.aperiodic.copy()
        self.gamma = self.aperiodic.copy()
        self.r_sq =  []
        self.fit_error = []


    def fitFOOOF(self,bg_param = 'knee',theta_range = (3,8),alpha_range = (8,12),beta_range = (12,30),
                      gamma_range = (30,40)):


        # clear results from any previous fit
        self.clearFOOOF()

        # Initialize FOOOF model
        fm = FOOOF(background_mode=bg_param) # using knee parameter leads to a better fit over long freq ranges

        # Fit FOOOF model
        fm.fit(self.freqs,self.psd_mean,self.frange)

        # access fit attributes
        #print('Background params [intercept, (knee), slope]:',fm.background_params_,'\n')
        #print('Peak params [CF: Center Frequencies, Amp: Amplitude Values, BW: Bandwidths]: \n \n',fm.peak_params_,
        #     '\n')
        #print('R Squared of full model fit to log10(psd):', fm.r_squared_,'\n')
        #print('root mean sq error of model fit', fm.error_,'\n')

        #collect model fit results (bg params, peak params, r sq, fit error, gaussian params)
        self.aperiodic,peak_params,self.r_sq,self.fit_error,gauss_params = fm.get_results()

        # parse peak_params into narrowbands
        # loop through peaks
        for i in (arange(0,peak_params.shape[0])):
            cf = peak_params[i][0]
            if (cf>=theta_range[0])&(cf<=theta_range[1]):
                self.theta = peak_params[i]
            elif (cf>alpha_range[0])&(cf<=alpha_range[1]):
                self.alpha = peak_params[i]
            elif (cf>beta_range[0])&(cf<=beta_range[1]):
                self.beta = peak_params[i]
            elif (cf>gamma_range[0])&(cf<=gamma_range[1]):
                self.gamma = peak_params[i]

        #collect fooof model
        self.fm = fm

    def plotFOOOF_fits(self,thisAx=None):
        if thisAx==None:
            fig = figure()
            thisAx = subplot(111)

        # plot fits
        self.fm.plot(ax=thisAx)

        # title is r-sq
        thisAx.set_title('r-sq = '+str(np.around(self.r_sq,decimals=2)),loc='right')

        # reset legend
        thisAx.legend()

    def plotFOOOF_params(self,ax=None,kol='blue'):
        if ax==None:
            fig = figure()
            ax = subplot(111)
        #set transparency
        alfa = 0.8

        # aperiodic slope
        plot(0,self.aperiodic[2],'o',color=kol,alpha = alfa)

        # for each narrowband (CF vs. amplitude)
        plot(self.theta[0],self.theta[1],'o',color=kol,alpha = alfa)
        plot(self.alpha[0],self.alpha[1],'o',color=kol,alpha = alfa)
        plot(self.beta[0],self.beta[1],'o',color=kol,alpha = alfa)
        plot(self.gamma[0],self.gamma[1],'o',color=kol,alpha = alfa)

        xlabel('Peak frequency')
        ylabel('Peak Amplitude')
        ax.set_xticks([0,3,8,12,30,40])
        ax.set_xticklabels(['SLOPE','3','8','12','30','40'])

    ### Time frequency analyses
    ### POWER
    def calcPow(self,buff_ms,frange,numfreqs,HFA_range = (70,120),theta_range = (3,8),alpha_range = (8,12),beta_range = (12,30),
                      gamma_range = (30,40)):
        """This function calculates log-transformed and normalized power surrounding the event_type.
        !!!!Note that normalization is only done within this specific event type so these power values
        cannot be used to compare power across trial types, but instead should only be used to relate to
        phase data related to these events

        Inputs:
        buff_ms: buffer associated with the EEG clip. This gets cropped out after wavelet is performed
        frange: frequency range e.g, (3,40)
        numfreqs: number of wavelets to fit within the freq range e.g., 35 fits 35 wavelets

        freq ranges for binning into power vectors
        HFA_range = (70,120),theta_range = (3,8),alpha_range = (8,12),beta_range = (12,30),
                          gamma_range = (30,40)



        Updates self with the following data:
        freqs_wave: frequency vector, useful for plotting
        zpowmat: log-transformed, and z-scored power matrix (see note about normalization routine above)
        zHFA:


        """
        # calc SR and freqs
        sr = self.eeg_clip.info['sfreq']
        myfreqs = np.logspace(np.log10(frange[0]), np.log10(frange[1]), num=numfreqs)

        # get power
        power = time_frequency.tfr_array_morlet(self.eeg_clip,sfreq=sr, freqs=myfreqs, n_cycles=5,output='power')

        # convert data from this object to a matrix, dropping the singleton dim, to get trials x freq x time
        powmat = np.squeeze(power.data)

        # drop the buffer
        powmat = powmat[:,:,int(buff_ms/1000.*sr):-1*int(buff_ms/1000.*sr)]

        # log transform the power
        powmat = log10(powmat)

        # calculate mean and std within freq, across time and trials
        mean_pow = powmat.mean(axis=0).mean(axis=1)
        std_pow = powmat.std(axis=0).std(axis=1)

        # reshape these arrays to have 3 dimensions (need to do this before next ste)
        mean_pow = np.reshape(mean_pow,(1,np.shape(powmat)[1],1))
        std_pow = np.reshape(std_pow,(1,np.shape(powmat)[1],1))

        # repeat these values across trials and across time
        mean_pow = np.repeat(np.repeat(mean_pow,np.shape(powmat)[0],axis=0),np.shape(powmat)[2],axis=2)
        std_pow = np.repeat(np.repeat(std_pow,np.shape(powmat)[0],axis=0),np.shape(powmat)[2],axis=2)

        #z-score power within-freq (across trials)
        z_powmat = (powmat-mean_pow)/std_pow

        # calculate time_window_samples (assoc w z_powmat after buffer removed)
        time_window_samp = linspace(int(self.time_window_ms[0]/1000*sr),int(self.time_window_ms[1]/1000*sr),np.shape(z_powmat)[2])

        # bin into specific ranges (results are trials x time)
        z_hfa =  np.mean(z_powmat[:,(myfreqs>=HFA_range[0])&(myfreqs<HFA_range[1]),:],axis=1).squeeze()
        z_theta =  np.mean(z_powmat[:,(myfreqs>=theta_range[0])&(myfreqs<theta_range[1]),:],axis=1).squeeze()
        z_alpha =  np.mean(z_powmat[:,(myfreqs>=alpha_range[0])&(myfreqs<alpha_range[1]),:],axis=1).squeeze()
        z_beta =  np.mean(z_powmat[:,(myfreqs>=beta_range[0])&(myfreqs<beta_range[1]),:],axis=1).squeeze()
        z_gamma =  np.mean(z_powmat[:,(myfreqs>=gamma_range[0])&(myfreqs<gamma_range[1]),:],axis=1).squeeze()

        # update self
        self.freqs_wave = myfreqs
        self.time_window_samp = time_window_samp
        self.HFA_range = HFA_range
        self.theta_range = theta_range
        self.beta_range = beta_range
        self.alpha_range = alpha_range
        self.gamma_range = gamma_range
        self.z_powmat = z_powmat
        self.z_HFA =  z_hfa
        self.z_theta =  z_theta
        self.z_alpha =  z_alpha
        self.z_beta =  z_beta
        self.z_gamma =  z_gamma

    def plotZPowMat(self,ax = None):
        """This function plots mean power. It assumes that self contains a power matrix that has already been
        log-transformed and z-scocred within frequency in the format ntrials x nfreq x time with buffers already cropped"""
        if ax==None:
            fig = figure()
            ax = subplot(111)
        fig = gcf()
        cax = ax.matshow(self.z_powmat.mean(axis=0).squeeze(), aspect='auto', cmap='RdBu_r');
        ax.invert_yaxis()
        yticks(np.arange(0, len(self.freqs_wave)), np.round(self.freqs_wave), fontsize=12)
        xlabel('Time (ms)', fontsize=14); ylabel('Frequency (Hz)', fontsize=14)
        cb = colorbar(cax); cb.set_label('log Power')
        ax.set_xticks(linspace(0,np.shape(self.z_powmat)[2],5))
        xvals_str = np.round(linspace(self.time_window_ms[0],self.time_window_ms[1],5)).astype('str')
        ax.set_xticklabels(xvals_str)

    def plotZPowVec(self,ax = None):
        C = cm.get_cmap('tab10')
        alfa = 0.65
        xval = linspace(self.time_window_ms[0],self.time_window_ms[1],np.shape(self.z_HFA)[1])
        hfa = plot(xval,self.z_HFA.mean(axis=0),color=C.colors[0],label='HFA'+str(self.HFA_range),alpha=alfa)
        theta = plot(xval,self.z_theta.mean(axis=0),color=C.colors[1],label='theta'+str(self.theta_range),alpha=alfa)
        beta = plot(xval,self.z_beta.mean(axis=0),color=C.colors[2],label='beta'+str(self.beta_range),alpha=alfa)
        vlines(0,ax.get_ylim()[0],ax.get_ylim()[1]);
        xlabel('Time from '+self.event_type+' (ms)')
        ax.legend()
        ax.tick_params(axis = 'x',direction = 'out', top = True, bottom = False, labeltop = True, labelbottom = False)
        ax.xaxis.set_label_position('top')


    ### PHASE
    def calcPhase(self,buff_ms,frange,numfreqs):
        """This function calculates phase. """

        # calc SR, and freq vector
        sr = self.eeg_clip.info['sfreq']
        myfreqs = np.logspace(np.log10(frange[0]), np.log10(frange[1]), num=numfreqs)

        # calc phase
        all_phase = time_frequency.tfr_array_morlet(self.eeg_clip,sfreq=sr, freqs=myfreqs, n_cycles=5,output='phase')

        # convert data from this object to a matrix, dropping the singleton dim, to get trials x freq x time
        all_phase = np.squeeze(all_phase)

        # drop the buffer
        all_phase = all_phase[:,:,int(buff_ms/1000.*sr):-1*int(buff_ms/1000.*sr)]

        # calculate time_window_samples (assoc w all_phase after buffer removed)
        time_window_samp = linspace(int(self.time_window_ms[0]/1000*sr),int(self.time_window_ms[1]/1000*sr),np.shape(all_phase)[2])


        # update self
        self.freqs_wave = myfreqs
        self.time_window_samp = time_window_samp
        self.all_phase = all_phase

    def calcPhaseDiff(self,other):
        """This function takes two instances of electrodes that have all_phase computed using
        the same parameters. It computes the phase difference between the phase distributions"""

        if np.shape(self.all_phase) == np.shape(other.all_phase):
            self.all_phase_diff =circ.cdiff(self.all_phase,other.all_phase)
            self.other_uElbl = other.uElbl
            self.other_roi = other.metadata['ROI'].iloc[0]
        else:
            print('all phase ')

    def calcPLV(self,freq_interest,time_win_interest,usePhaseDiff=False):
        """
        This function calculates the PLV across trials for a specified frequency and time segment
        inputs:
        freq_interest = (3,3)
        time_win_interest = (-500,-300)

        Required the following in self
        self.all_phase = ndarray of phase values (trials x freqency x time (samples))
        self.freqs_wave = array of frequencies used to calculate phase distributoin
        self.time_windows_samp = samples associated with phase array (relataive to cue onset)

        # if usePhaseDiff == True
        it will use the phase distribution stored in self.phase_diff which represents the phase difference
        between self and other (another instance of Electrode)

        """

        # calculate sr
        sr = self.eeg_clip.info['sfreq']

        #
        if usePhaseDiff == True:
            all_phase = self.all_phase_diff
        else:
            all_phase = self.all_phase

        #freq_interest = (3,3)
        f_start = np.argmin(np.absolute(self.freqs_wave-freq_interest[0]))
        f_end = np.argmin(np.absolute(self.freqs_wave-freq_interest[1]))

        #time_win_interest = (-500,-300)
        samp_start = np.argmin(np.absolute(self.time_window_samp-int((time_win_interest[0]/1000*sr))))
        samp_end = np.argmin(np.absolute(self.time_window_samp-int((time_win_interest[1]/1000*sr))))


        # get phase segment and do circ mean
        if f_start==f_end: # only one frequency
            phs = all_phase[:,f_start,samp_start:samp_end]
            # do mean across time range and freq range
            phs_avg = circ.mean(phs,axis=1)
        else:
            phs = all_phase[:,f_start:f_end,samp_start:samp_end]
                # do mean across time range and freq range
            phs_avg = circ.mean(circ.mean(phs,axis=2),axis=1)

        # calc resultant vector length anad Raleigh statistic
        rvl = circ.resultant_vector_length(phs_avg)
        p,z = circ.rayleigh(phs_avg);

        #update self
        self.PLV_freq_interest = freq_interest
        self.PLV_f_start = f_start
        self.PLV_f_end = f_end
        self.PLV_time_win_interest = time_win_interest
        self.PLV_p = p
        self.PLV_z = z
        self.PLV_rvl= rvl
        self.phs_avg = phs_avg


    def calcPLV_over_time(self, freq_interest,twin=None,tstep=200,usePhaseDiff=False):
        """ This function computes phase locking value over non-overlapping time windows for a given phase distribution.
        Can be a distribution of phases for single-electrode PLV or phase-differences for pairwise PLV
        Inputs:

        freq_interest = (3,8) # frequency range of interest
        twin = (-1800,1800) # time range over which to loop over, must be within the range of phase distribution
        tstep = 200 # length of time window

        # self must contain the following:
        self.all_phase # phase distribution
        self.freqs_wave # frequency vector used for phase distribution
        self.time_windows_samp = samples associated with phase array (relataive to cue onset)

        # if usePhaseDiff == True
        it will use the phase distribution stored in self.phase_diff which represents the phase difference
        between self and other (another instance of Electrode)

        p_thresh = 0.05 # p-theshold to highlight

        """
        #calc sr
        sr = self.eeg_clip.info['sfreq']

        # calc default time windows
        if twin == None:
            twin = self.time_window_ms


        # calculate start times for each time window
        t_starts = arange(twin[0],twin[1],tstep)

        #Initialize containers
        z_vec = np.empty(len(t_starts));
        p_vec = np.empty(len(t_starts));
        t_vec = np.empty(len(t_starts));

        count = -1 # count
        for t in t_starts:
            count+=1
            time_win_interest = (t,t+tstep);#(-1000,-800)
            t_vec[count] = np.mean(time_win_interest)

            # run PLV on this time window
            self.calcPLV(freq_interest,time_win_interest,usePhaseDiff);

            # collect stats
            p_vec[count] = self.PLV_p
            z_vec[count] = self.PLV_z

        #update self
        self.PLVot_t_vec = t_vec
        self.PLVot_p_vec = p_vec
        self.PLVot_z_vec = z_vec
        self.PLVot_freq_interest = freq_interest
        self.PLVot_time_win_interest = twin
        self.PLVot_tstep = tstep
        # store the time window of max PLV (for polar plotting)
        self.PLVot_maxZ_time_win = (t_starts[np.argmax(z_vec)],t_starts[np.argmax(z_vec)]+tstep)


    def plotPLV_polar(self,ax = None):
        """This function plots a polar histram of currently averaged phase distribution"""
        if ax==None:
            fig = figure()
            ax = subplot(111)

        ax.polar=True
        hist(self.phs_avg, 40);
        title('time '+str(self.PLV_time_win_interest)+' z = '+str(np.around(self.PLV_z,decimals=2))+' p = '+str(np.around(self.PLV_p,decimals=2)))

    def plotPLV_over_time(self,ax = None,p_thresh=0.05):
        """This function plots PLV results over time"""

        if ax==None:
            fig = figure()
            ax = subplot(111)

        plot(self.PLVot_t_vec,self.PLVot_z_vec)
        plot(self.PLVot_t_vec[self.PLVot_p_vec<=p_thresh],self.PLVot_z_vec[self.PLVot_p_vec<=p_thresh],'or')
        xlabel('Time from cue onset (ms)')
        ylabel('Raleigh (z)')
        title('Phase locking value for Frequency '+str(self.PLVot_freq_interest))
