#import packages
import numpy as np # numpy
import pandas as pd # pandas
import mne
import pickle
import os
from scipy.io import loadmat # to load matlab
from scipy import stats,ndimage,signal,spatial
import numpy.random as rand
import pycircstat as circ
import fooof as ff
import tensorpac as tp

# plotting
import matplotlib
#matplotlib.use('macosx')
matplotlib.use('Qt5Agg')

# this makes tight axes 
matplotlib.rcParams['axes.autolimit_mode'] = 'round_numbers'
matplotlib.rcParams['axes.xmargin'] = 0
matplotlib.rcParams['axes.ymargin'] = 0
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt 

# load dataset. This will return a list of dictionaries
fpath = '/Users/ashwinramayya/Dropbox/neurosurgery_research/data/scratch/clusStruct/'
fname = 'groupDataStruct.mat'
matFile =loadmat(fpath+fname,squeeze_me=True)

# convert to a dataframe. This can be accessed by all functions below
groupDf = pd.DataFrame(matFile['groupDataStruct'])

# force numeric datatypes (for plotting purposes); do not include arrays
ncol_list = ['clusId', 'sr', 'postInst_HFA', 'postMove_HFA', 'effectorSelective','instructionSelective', 'x', 'y', 'z']

for n in ncol_list:
          groupDf[n] = pd.to_numeric(groupDf[n],errors = 'raise')


# create scratch dir
scratch_dir = os.getcwd()+'/scratch/'
if os.path.exists(scratch_dir)==False:
     os.mkdir(scratch_dir)

# create fig dir
fig_dir = os.getcwd()+'/figs/'
if os.path.exists(fig_dir)==False:
     os.mkdir(self.fig_dir)


##############################################################################
#ELECTRODE CLASS
###############################################################################

class Electrode:
          # performs analyses on data from a single electrode (across all sessions), include pairwise analyses between (self,other)

          # constructor
          def __init__(self,uElbl,paramsDict = None):
              #Input
              # subj Id
              # paramsDict = Optional. dictionary of params to update. Will overwrite default param values

              # get indx
              self.e_idx = groupDf.eval('uElbl==@uElbl').to_numpy().nonzero()[0][0]          
              # get dictionary representing this electrode's data
              # We will continue to update this series with results from computations below
              self.elec_dict = groupDf.loc[self.e_idx].to_dict()

              # params
              self.params = {'tmin_ms':-1*self.elec_dict['eegConfig']['priorMS'], 'tmax_ms':self.elec_dict['eegConfig']['postMS'], # 
                                 'eeg_notch_filter':False,
                                 'eeg_filt_l_freq':62,
                                 'eeg_filt_h_freq':58,
                                 'buffer_ms':self.elec_dict['eegConfig']['bufferMS'],
                                 'frange_HFA':(70,125),
                                 'frange_LFA':(2,50),
                                 'frange_thetaL':(3,6),
                                 'frange_thetaH':(6,8),
                                 'frange_thetaW':(2,10),
                                 'frange_theta':(3,8),
                                 'frange_alpha':(8,12),
                                 'frange_beta':(12,30),
                                 'frange_gamma':(30,50),
                                 'wave_frange':(3,200),
                                 'wave_number_of_freq':6,
                                 'trange_periEvent':500,
                                 'trange_periEventShort':250}

              # add wave number of cycles
              self.params['wave_number_of_cycles'] = 5
              #self.params['wave_number_of_cycles'] = np.flip(np.linspace(3,7,self.params['wave_number_of_freq']))

              # update this dict if dictionary values are provided
              if paramsDict != None:
                         self.params.update(paramsDict)

          ########################## 
          #GENERAL FUNCTIONS
          ##########################
          
          def updateParamsDict(self,key,val):
              # overwrites similar function in Subject class
              # update SE params dictionary
              self.params[key] = val

          def initResultsDict(self,other = None,incl_LFO = False):

                # This function returns a dictionary populated with key data from elec_dict

                res_dict = {}    


                # elec_dict data for convenience
                res_dict['subj'] = self.elec_dict['subj']
                res_dict['uElbl'] = self.elec_dict['uElbl']
                res_dict['clusId'] = self.elec_dict['clusId']
                res_dict['ROI'] = self.elec_dict['ROI']
                res_dict['postInst_HFA'] = self.elec_dict['postInst_HFA']
                res_dict['postMove_HFA'] = self.elec_dict['postMove_HFA']
                res_dict['postInst_theta'] = self.elec_dict['postInst_theta']
                res_dict['postMove_theta'] = self.elec_dict['postMove_theta']
                res_dict['postInst_beta'] = self.elec_dict['postInst_beta']
                res_dict['postMove_beta'] = self.elec_dict['postMove_beta']
                res_dict['postInst_LFA'] = self.elec_dict['postInst_LFA']
                res_dict['postMove_LFA'] = self.elec_dict['postMove_LFA']
                res_dict['effectorSelective'] = self.elec_dict['effectorSelective']
                res_dict['instructionSelective'] = self.elec_dict['instructionSelective']
                res_dict['x'] = self.elec_dict['x']
                res_dict['y'] = self.elec_dict['y']
                res_dict['z'] = self.elec_dict['z']

                if incl_LFO == True:
                     self.calc_lfo(ev_type=None,offset_ms=0, duration_ms = 5000)
                     res_dict['lfo_cf'] = self.lfo_dict['LFO'][0]
                     res_dict['theta_cf'] = self.lfo_dict['theta'][0]
                     res_dict['beta_cf'] = self.lfo_dict['beta'][0]




                # if there is an other object
                if other != None:
                    res_dict['other_subj'] = other.elec_dict['subj']
                    res_dict['other_uElbl'] = other.elec_dict['uElbl']
                    res_dict['other_clusId'] = other.elec_dict['clusId']
                    res_dict['other_ROI'] = other.elec_dict['ROI']
                    res_dict['other_postInst_HFA'] = other.elec_dict['postInst_HFA']
                    res_dict['other_postMove_HFA'] = other.elec_dict['postMove_HFA']
                    res_dict['other_effectorSelective'] = other.elec_dict['effectorSelective']  
                    res_dict['other_instructionSelective'] = other.elec_dict['instructionSelective']
                    res_dict['x'] = other.elec_dict['x']
                    res_dict['y'] = other.elec_dict['y']
                    res_dict['z'] = other.elec_dict['z']

                    # pairwise anatomical distance 
                    xyz = [[self.elec_dict['x'],self.elec_dict['y'],self.elec_dict['z']],[other.elec_dict['x'],other.elec_dict['y'],other.elec_dict['z']]]
                    res_dict['pairwise_dist'] = spatial.distance.pdist(xyz)

                    if incl_LFO == True:
                         other.calc_lfo(ev_type=None,offset_ms=0, duration_ms = 5000)
                         res_dict['other_lfo_cf'] = other.lfo_dict['LFO'][0]
                         res_dict['other_theta_cf'] = other.lfo_dict['theta'][0]
                         res_dict['other_beta_cf'] = other.lfo_dict['beta'][0]




                return res_dict

          def save_pickle(self,obj, fpath, verbose=True):
              """Save object as a pickle file. Written by Dan Schonhaut. 11/2018"""
              with open(fpath, 'wb') as f:
                         pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
              if verbose:
                         print('Saved {}'.format(fpath))

          def load_pickle(self,fpath):
              """Return object. Written by Dan Schonhaut. 11/2018"""
              with open(fpath, 'rb') as f:
                         obj = pickle.load(f)
              return obj

          # drop buffer
          def dropBuffer(self,datamat_2d):
              # This is a general function to drop the buffer from a 2d data matrix.
              # Inputs:
              # datamat_2d ... any 2d data matrix (rows = trials, columns = time). Assumes the data are in samples and range from (-bufferMS+tmin_ms):(bufferMS+tmax_ms). can provide ERP, plv_xtrial_2d as inputs

              # get buffer index
              buff_idx = int((self.params['buffer_ms']/1000)*self.elec_dict['sr'])

              # if 1 dim array is passed in, add an axis 
              if datamat_2d.ndim==1:
                         datamat_2d = np.array(datamat_2d[np.newaxis,:])
                         datamat_2d_nobuff = datamat_2d[:,buff_idx:-buff_idx];
                         datamat_2d_nobuff = np.squeeze(datamat_2d_nobuff) # back to orignial dim
              else:
                         datamat_2d_nobuff = datamat_2d[:,buff_idx:-buff_idx];

              return datamat_2d_nobuff

          ## ms_to_samples
          def ms_to_samples(self,ms):
              # Input
              #ms ... ms value to convert to samples (int or array)
              samples = (ms/1000)*self.elec_dict['sr']
              return samples

          # samples_to_ms
          def samples_to_ms(self,samples):
              # Input
              #samples ... sample value to convert to samples (int or array)
              ms = (samples/self.elec_dict['sr'])*1000
              return ms


          def makeTimeBins(self,time_bin_size_ms = 100):
              # This function generates a time window dictionary and updates self object. Can be used for additional computations

              #calc time bins
              time_bin_starts = np.arange(self.params['tmin_ms'],self.params['tmax_ms'],time_bin_size_ms)

              #create time bin dictionary
              self.time_bin_dict = {'size_ms':time_bin_size_ms,
                                                                  'starts_ms':time_bin_starts,
                                                                  'lbl':(time_bin_starts+time_bin_size_ms/2).astype('int')}

          def clip_data(self,datamat_nobuff, offset_ms = 0, duration_ms = 1000 ,xval_ms = None):
                # This function clips data based on a time window. 
                # Inputs:
                # datamat_nobuff ... a matrix with time series data. Can be array like (1d) or a 2d matrix (rows = trials, col = time). It expects data without a buffer
                # xval_ms ... array indicating time stamps relative to event for the data matrix (if None, it creates this array based on tmin and max in params data)
                # offset_ms ... clip start relative to event onset 
                # duration_ms ... duration of clip in ms

                # if no xval_ms is given,assume we are dealing with the entire eeg clip
                if xval_ms == None:
                     xval_ms = np.linspace(self.params['tmin_ms'],self.params['tmax_ms'],np.shape(datamat_nobuff)[-1])

                # find samp start and stop
                samp_start_idx = np.argmin(np.abs(xval_ms - offset_ms))
                samp_end_idx = np.argmin(np.abs(xval_ms - (offset_ms+duration_ms)))

                xval_ms_clip = xval_ms[samp_start_idx:samp_end_idx]

                # check if we are dealign with a 1 d or 2d array
                if datamat_nobuff.ndim == 1:
                     datamat_clip = datamat_nobuff[samp_start_idx:samp_end_idx]
                elif datamat_nobuff.ndim ==2:
                     datamat_clip = datamat_nobuff[:,samp_start_idx:samp_end_idx]

                return datamat_clip, xval_ms_clip


          def data2timeBins(self,datamat_2d, time_bin_size_ms = 100, isCircData = False):
                     # This function takes a 2D power matrix (trials x samples) and bins it into time windows
                     # returns the binned matrix as a variable (does not store in self)
                     # Inputs
                     # datamat_2d..2D matrix (trials x samples)
                     # time_bin_size_ms = 100....length of time window to bin
                     # isCircData .. if True, it computes circular mean for each time bin

                     # if 1 dim array is passed in, add an axis 
                     if datamat_2d.ndim==1:
                                datamat_2d = np.array(datamat_2d[np.newaxis,:])
                                squeeze_at_end = True
                     else:
              	          squeeze_at_end = False


                     # make time bins
                     self.makeTimeBins(time_bin_size_ms=time_bin_size_ms)

                     # get sample vector
                     samps_vec = np.arange(int(self.ms_to_samples(self.params['tmin_ms'])),
                                                                  (self.ms_to_samples(self.params['tmax_ms'])))

                     # initialize bin vector
                     datamat_2d_ms = np.empty((datamat_2d.shape[0],self.time_bin_dict['starts_ms'].shape[0]))

                     # loop through time bins
                     for t in np.arange(0,len(self.time_bin_dict['starts_ms'])):
                                samp_start = int((self.time_bin_dict['starts_ms'][t]/1000)*self.elec_dict['sr'])
                                samp_start_idx = np.argmin(np.abs(samps_vec-samp_start))

                                samp_end = int(samp_start+(self.time_bin_dict['size_ms']/1000)*self.elec_dict['sr'])
                                samp_end_idx = np.argmin(np.abs(samps_vec-samp_end))

                                if isCircData == True:
                                          datamat_2d_ms[:,t] = circ.mean(datamat_2d[:,samp_start_idx: samp_end_idx],axis=1)
                                else:
                                          datamat_2d_ms[:,t] = np.nanmean(datamat_2d[:,samp_start_idx:samp_end_idx],axis=1)

                     # return to original dimensions
                     if squeeze_at_end == True:
                                datamat_2d_ms = datamat_2d_ms.squeeze()

                     return datamat_2d_ms

          # get continous timeseries of EEG across trials
          def catEEG(self, ev_type = None, ev_item = None,offset_ms = 0, duration_ms = 5000):
                     #This function concatonates eeg data across trials to generate a  single vector. It adds a buffer at the beginning and the end that should be dropped after running the subsequent analysis using (self.dropBuffer())
                     # Inputs:
                     # ev_type....if = None,it will use all event type (wait, move,  instruct); otherwise it will only use data from the specified event type
                     # ev_item...option to filter by ev_item
                     # offset_ms...sets the start of the clip for each trial (0 starts the  clip from cue onset, -1000 ms starts 1 sec prior, 1000 would start 1 sec after stim onset)
                     # duration_ms ...duration length of the clip to take from each session (5000 means clip 5 sec from offset_ms)

                     # returns 
                     # eeg_vec_w_buff .... concatenated eeg data (single time series) with buffer at start and end

                     # get eeg matrix
                     eeg = self.elec_dict['subjEEG']

                     # get buffer from beginning of first trial and last trial; will add this to the concatenated time series
                     buff_idx = int((self.params['buffer_ms']/1000)*self.elec_dict['sr'])
                     buff_start = eeg[0,:buff_idx]
                     buff_end = eeg[-1,-buff_idx:]
 
                     # drop the buffer (so we can calculate xrange of clips)
                     eeg_noBuff = self.dropBuffer(eeg)


                     # find clip start and stop indices (for each trial)
                     samp_xval = np.linspace(self.ms_to_samples(self.params['tmin_ms']),self.ms_to_samples(self.params['tmax_ms']),np.shape(eeg_noBuff)[1])
                     samp_start_idx = np.argmin(np.abs(samp_xval-self.ms_to_samples(offset_ms)))
                     samp_end_idx = np.argmin(np.abs(samp_xval-self.ms_to_samples(offset_ms+duration_ms)))

                     # filter eeg
                     # by event type
                     if ev_type == None:
                                type_idx = np.ones((len(self.elec_dict['ev_type']))).astype('bool')
                     else:
                                type_idx = self.elec_dict['ev_type'] == ev_type

                     # by event item
                     if ev_item == None:
                                item_idx = np.ones((len(self.elec_dict['ev_item']))).astype('bool')
                     else:
                                item_idx = self.elec_dict['ev_item'] == ev_item


                     # filtered eeg by trial and sample (according to offset_ms and  duration_ms)
                     eeg_noBuff_filt = eeg_noBuff[type_idx&item_idx,samp_start_idx:samp_end_idx]


                     # reshape to a 1 d array by inferring length from the matrix
                     eeg_vec = np.reshape(eeg_noBuff_filt, newshape = -1)

                     # add buffer
                     eeg_vec_w_buff = np.concatenate((buff_start,eeg_vec,buff_end)) 

                     # return 
                     return eeg_vec_w_buff

          def eeg2hilb(self,eeg_w_buff,fRange):
                #This function filters within fRange and applies the hilbert transfrom on eeg data. It expects buffered data. It drops the buffer before returning phase and power data
                # Inputs: 
                # eeg_w_buff canvstack of timeseries organized by trial; or a single time series. NOTE!!!:: It runs ~100 times faster (.15 s v.s. 15 s) if it is given a matrix rather than a single concatonated time series
                # NOTE: Returns phase in RADIANS

                # convert to float
                eeg_w_buff = eeg_w_buff.astype('float')

                # filter data
                eeg_filt_w_buff = mne.filter.filter_data(data = eeg_w_buff, sfreq = self.elec_dict['sr'],l_freq = fRange[0],h_freq = fRange[1],method='fir',verbose=False)

                # apply hilbert, returns a complex valued vector
                eeg_hilbert_w_buff = signal.hilbert(eeg_filt_w_buff)


                # drop the buffer
                eeg_hilbert = self.dropBuffer(eeg_hilbert_w_buff)


                # compute phase
                hilb_phase = np.angle(eeg_hilbert,deg=False)

                # compute power
                hilb_power = np.abs(eeg_hilbert)**2

                # log transform power
                hilb_power = np.log10(hilb_power)

                return hilb_phase,hilb_power

          def eeg2wave(self,eeg_w_buff,fRange,numfreqs=None,numcycles=None):
                #This function filters within fRange and applies the hilbert transfrom on eeg data. It expects buffered data. It drops the buffer before returning phase and power data
                # Inputs: 
                # eeg_w_buff can be vstack of timeseries organized by trial; or a single time series. 
                # fRange .. freq range for wavelets
                # numfreqs ... number of freqs (num wavelets)
                # numcycles .. num of cycles (~wave number)

                # parse inputs (if none, it will default to params)
                if numfreqs==None:
                     numfreqs = self.params['wave_number_of_freq']
                if numcycles == None:
                     numcycles = self.params['wave_number_of_cycles']

                myfreqs = np.logspace(np.log10(fRange[0]), np.log10(fRange[1]),num=numfreqs)


                # convert to float
                eeg_w_buff = eeg_w_buff.astype('float')

                # reshape into (n_epochs,n_chans,n_times)
                eeg_w_buff = eeg_w_buff[:,np.newaxis,:]

                # get complex data
                wave_complex = mne.time_frequency.tfr_array_morlet(eeg_w_buff,sfreq=self.elec_dict['sr'].astype('float'), freqs=myfreqs, n_cycles=numcycles,output='complex')

                # reshape to freq x time x trials
                wave_complex = np.moveaxis(wave_complex,(2,3,0),(0,1,2)).squeeze()

                # drop the buffer (manually bc this is 3D data)
                wave_complex = wave_complex[:,int(self.ms_to_samples(self.params['buffer_ms'])):-1*int(self.ms_to_samples(self.params['buffer_ms'])),:]
 
                # compute phase
                wave_phase = np.angle(wave_complex)
 
                # calculate power
                wave_power = (wave_complex * wave_complex.conj()).real

                # log transform power
                wave_logpower = np.log10(wave_power)

                return wave_phase,wave_logpower,wave_complex
          def phsMat2timeBins(self,phsMat_no_buff, time_bin_size_ms = 100):
              # This function takes a 2D phase matrix (trials x samples) and bins it into time windows
              # returns the binned matrix as a variable (does not store in self)
              # Inputs
              # phsMat..2D matrix (trials x samples).
              # time_bin_size_ms = 100....length of time window to bin

              # make time bins
              self.makeTimeBins(time_bin_size_ms=time_bin_size_ms)

              # get sample vector
              samps_vec = np.arange(int(self.ms_to_samples(self.params['tmin_ms'])),
                                                                  (self.ms_to_samples(self.params['tmax_ms'])))

              # initialize bin vector
              phsMat_bin_no_buff = np.empty((phsMat_no_buff.shape[0],self.time_bin_dict['starts_ms'].shape[0]))

              # loop through time bins
              for t in np.arange(0,len(self.time_bin_dict['starts_ms'])):
                         samp_start = self.ms_to_samples(self.time_bin_dict['starts_ms'][t])
                         samp_start_idx = np.argmin(np.abs(samps_vec-samp_start))

                         samp_end = int(samp_start+self.ms_to_samples(self.time_bin_dict['size_ms']))
                         samp_end_idx = np.argmin(np.abs(samps_vec-samp_end))

                         phsMat_bin_no_buff[:,t] = circ.mean(phsMat_no_buff[:,samp_start_idx:samp_end_idx],axis=1)

              return phsMat_bin_no_buff

          def calc_phase(self,other = None, fRange = (3,8), method = 'filt_hilb',wave_numfreqs=None,wave_numcycles=None):
                     # This function calculates phase data for the electrode for for all trials for the given freq range using the filter/hilbert method (it has a method in place to use buffered data to calculate phase and drops the buffer after). If other is another Electrode object, it computes pairwise phase differences. It populates self with phsmat (intended to be used as a common computation in other analyses) 
                     # get eeg array
                     eeg_self = self.elec_dict['subjEEG']
                         
                     # calc phase data (drops buffer when returning phase and power data)
                     if method == 'filt_hilb':
                          phsmat_self, logpow_self = self.eeg2hilb(eeg_self,fRange)
                     elif method == 'wave':
                          # this returns a 3d mat (freq x time x shape). 
                          phsmat_self, logpow_self,wavecomp = self.eeg2wave(eeg_self,fRange = fRange,numfreqs=wave_numfreqs,numcycles=wave_numcycles)

                          #Collapse to 2d matrix (trials x time)
                          phsmat_self = np.swapaxes(circ.mean(phsmat_self,axis = 0),0,1)

                     # parse "other"
                     if other != None:
                                #if other is another electrode object. phsmat is a matrix of phase difference

                                if (other.elec_dict['subj']!=self.elec_dict['subj'])|(other.elec_dict['sr']!=self.elec_dict['sr']):
                                          raise NameError('OTHER must be an electrode from the same subject')

                                # calculate eeg from other (should have the same sample rate)
                                eeg_other = other.elec_dict['subjEEG']

 
                                #calculate phase mat from other
                                if method == 'filt_hilb':
                                     phsmat_other, logpow_other = other.eeg2hilb(eeg_other,fRange)
                                elif method == 'wave':
                                     phsmat_other, logpow_other,wavecomp = other.eeg2wave(eeg_other,fRange = fRange,numfreqs=wave_numfreqs,numcycles=wave_numcycles)
                                     
                                     #Collapse to 2d matrix (trials x time)
                                     phsmat_other = np.swapaxes(circ.mean(phsmat_other,axis = 0),0,1)

                                #check if there is a mismatch in sessions between self and other
                                sess_list_self = np.unique(self.elec_dict['ev_sess'])
                                sess_list_other = np.unique(other.elec_dict['ev_sess'])
                                if len(sess_list_self) != len(sess_list_other):
                                     # filter phase mat to get matching trials
                                     sess_list_all = np.unique(np.concatenate((sess_list_self,sess_list_other)))
                                     idx_to_keep_self = np.zeros(len(self.elec_dict['ev_sess'])).astype('bool')
                                     idx_to_keep_other = np.zeros(len(other.elec_dict['ev_sess'])).astype('bool')

                                     for s in sess_list_all:
                                          if (s in sess_list_self)&(s in sess_list_other):
                                                idx_to_keep_self[self.elec_dict['ev_sess']==s] = True
                                                idx_to_keep_other[other.elec_dict['ev_sess']==s] = True

                                     # filter data
                                     # 2d matrix
                                     phsmat_self = phsmat_self[idx_to_keep_self,:]
                                     phsmat_other = phsmat_other[idx_to_keep_other,:]
                                     
                                     # other data
                                     phsmat_ev_type = self.elec_dict['ev_type'][idx_to_keep_self]
                                     phsmat_ev_item = self.elec_dict['ev_item'][idx_to_keep_self]
                                     phsmat_eeg_self = self.elec_dict['subjEEG'][idx_to_keep_self,:]
                                     phsmat_eeg_other = other.elec_dict['subjEEG'][idx_to_keep_other,:]

                                else:
                                     phsmat_ev_type = self.elec_dict['ev_type']
                                     phsmat_ev_item = self.elec_dict['ev_item']
                                     phsmat_eeg_self = self.elec_dict['subjEEG']
                                     phsmat_eeg_other = other.elec_dict['subjEEG']
 
                                # calculate phase differences between phase time series
                                phsmat =circ.cdiff(phsmat_self,phsmat_other)

                     else:
                                phsmat = phsmat_self
                                phsmat_ev_type = self.elec_dict['ev_type']
                                phsmat_ev_item = self.elec_dict['ev_item']
                                phsmat_eeg_self = self.elec_dict['subjEEG']
                                phsmat_eeg_other = None

                     
                     # populate phasemat dict 
                     self.phsmat_dict = self.initResultsDict(other = other)
                     self.phsmat_dict['phsmat'] = phsmat
                     self.phsmat_dict['fRange'] = fRange
                     self.phsmat_dict['method'] = method
                     self.phsmat_dict['wave_numfreqs'] = wave_numfreqs
                     self.phsmat_dict['wave_numcycles'] = wave_numcycles
                     self.phsmat_dict['ev_type'] = phsmat_ev_type
                     self.phsmat_dict['ev_item'] = phsmat_ev_item
                     self.phsmat_dict['eeg_self'] = phsmat_eeg_self
                     self.phsmat_dict['eeg_other'] = phsmat_eeg_other

                     # label whethere phsmat is a difference measure or raw phase
                     if other!=None:
                                self.phsmat_dict['isDiff'] = True
                                self.phsmat_dict['other_uElbl']= other.elec_dict['uElbl']
                                self.phsmat_dict['phs_self'] = phsmat_self
                                self.phsmat_dict['phs_other'] = phsmat_other
                     else:
                                self.phsmat_dict['isDiff'] = False
                                self.phsmat_dict['other_uElbl']= None
                                self.phsmat_dict['phs_self'] = phsmat_self
                                self.phsmat_dict['phs_other'] = np.empty(np.shape(phsmat_self))
                                self.phsmat_dict['phs_other'][:] = np.nan 
          def get_phase_filteredByTrial(self,other=None,fRange=(3,8), ev_type=None, ev_item=None, method = 'filt_hilb',wave_numfreqs = None, wave_numcycles = None):

                     # check whether we already have the appropriate phase data cached
                     calcPhsFlag = True
                     if hasattr(self,'phsmat_dict') == True:
                                # we have some phase data
                                # now check if it is appropriate
                                if (other == None) & (self.phsmat_dict['isDiff'] == False) & np.all((self.phsmat_dict['fRange']==fRange)) & (self.phsmat_dict['method']==method)& (self.phsmat_dict['wave_numfreqs']==wave_numfreqs)& (self.phsmat_dict['wave_numcycles']==wave_numcycles):
                                  calcPhsFlag = False
                                elif (other != None) & (self.phsmat_dict['other_uElbl'] == other.elec_dict['uElbl']) & np.all((self.phsmat_dict['fRange']==fRange)) & (self.phsmat_dict['method']==method)& (self.phsmat_dict['wave_numfreqs']==wave_numfreqs)& (self.phsmat_dict['wave_numcycles']==wave_numcycles):
                                  calcPhsFlag = False

                     # if we havent calculated the appropriate phase data, recalculate it here
                     if calcPhsFlag==True:
                          self.calc_phase(other = other,fRange = fRange,method = method,wave_numfreqs = wave_numfreqs, wave_numcycles = wave_numcycles)


                     # create a copy of phsmat 
                     phsmat_trFilt = np.copy(self.phsmat_dict['phsmat'])

                     # filter by event type
                     if ev_type == None:
                          type_idx = np.ones((len(self.phsmat_dict['ev_type']))).astype('bool')
                     else:
                          type_idx = self.phsmat_dict['ev_type'] == ev_type

                     # filter by event item
                     if ev_item == None:
                          item_idx = np.ones((len(self.phsmat_dict['ev_item']))).astype('bool')
                     else:
                          item_idx = self.phsmat_dict['ev_item'] == ev_item

                     # filter phsmat
                     phsmat_trFilt = phsmat_trFilt[type_idx&item_idx,:]             

                     return phsmat_trFilt

          def calc_null_phase(self, n_iters = 10000):
                # This function calculates a null phase distribution based on the cached phase distributions  a given electrode pair and time freq parameters. The goal is to maintain autocorrelation but remove pairwise correlations. For each iteration, it randomly selects a trial from phsmat_self and phsmat_other (with replacement) and circularly shifts by a random value. it computes the phase difference on null data and returns it

                # first check if we have computed true phase data 
                if hasattr(self,'phsmat_dict') != True:
                     raise NameError('calculate phase first')


                # then check if we already calculated null phase data
                calcNull = []
                if hasattr(self,'nullphsmat_dict') == True:        

                     # if we have null data, check if the method parameters match those of the true phase data
                     if np.all((self.nullphsmat_dict['fRange']==self.phsmat_dict['fRange'])) & (self.nullphsmat_dict['method']==self.phsmat_dict['method']) & (self.nullphsmat_dict['wave_numfreqs']==self.phsmat_dict['wave_numfreqs']) & (self.nullphsmat_dict['wave_numcycles']==self.phsmat_dict['wave_numcycles']) & (self.nullphsmat_dict['other_uElbl']==self.phsmat_dict['other_uElbl']):

                          # if so, return cached null phs data
                          calcNull = False
                          return self.nullphsmat_dict['null_phsmat']
                     else:
                          calcNull = True
                else:
                     calcNull = True
                if calcNull==True:
                     # compute null phase data
                     # load self and other
                     phs_self = self.phsmat_dict['phs_self']
                     phs_other = self.phsmat_dict['phs_other']

                     # init 
                     null_phs_self  = np.empty((n_iters,np.shape(phs_self)[1]))
                     null_phs_self[:] = np.nan
                     null_phs_other = np.empty((n_iters,np.shape(phs_other)[1]))
                     null_phs_other[:] = np.nan

                     # loop through iterations
                     for i in np.arange(0,n_iters):
                          

                          #randomly generate trial idx self and other
                          tr_idx_s = rand.randint(low=0, high=np.shape(phs_self)[0])
                          tr_idx_o = rand.randint(low=0, high=np.shape(phs_other)[0])

                          #randomly generate circ shift (roll) idx self and other
                          sh_idx_s = rand.randint(low=0, high=np.shape(phs_self)[1])
                          sh_idx_o = rand.randint(low=0, high=np.shape(phs_other)[1])

                          # populate null matrices by randomly circ shifting each trial
                          null_phs_self[i,:] = np.roll(phs_self[tr_idx_s,:],sh_idx_s)
                          null_phs_other[i,:] = np.roll(phs_other[tr_idx_o,:],sh_idx_o)


                          #print(i,tr_idx_s,tr_idx_o,sh_idx_s,sh_idx_o)


                     # calc phase difference
                     null_phsmat = circ.cdiff(null_phs_self,null_phs_other)


                     # save in dict for easy retrieval
                     self.nullphsmat_dict = {}
                     self.nullphsmat_dict['null_phsmat'] = null_phsmat
                     self.nullphsmat_dict['fRange'] = np.copy(self.phsmat_dict['fRange'])
                     self.nullphsmat_dict['method'] = np.copy(self.phsmat_dict['method'])
                     self.nullphsmat_dict['wave_numfreqs'] = np.copy(self.phsmat_dict['wave_numfreqs'])
                     self.nullphsmat_dict['wave_numcycles'] = np.copy(self.phsmat_dict['wave_numcycles'])
                     self.nullphsmat_dict['isDiff'] = np.copy(self.phsmat_dict['isDiff'])
                     self.nullphsmat_dict['other_uElbl']= np.copy(self.phsmat_dict['other_uElbl'])

                     return self.nullphsmat_dict['null_phsmat']
          def plot_scatter(self,x,y,ax = None,plotLine=True,polyfit=False):
              # general plotting function to plot a scatter plot for two variables and fit a line. Returns x and y
              #remove nans
              rm_bool = np.isnan(x) | np.isnan(y)
              x=x[rm_bool==False]
              y = y[rm_bool==False]
              # x vs. y

              if ax == None:
                f = plt.figure()
                ax = plt.subplot(111)
              plt.scatter(x,y,c='0.5',edgecolor = 'k',alpha =0.5)

              if plotLine==True:
                slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
                if polyfit == True:
                    polymod= np.poly1d(np.polyfit(x, y, 3))
                    x_mod = np.linspace(np.min(x),np.max(x),20)
                    plt.plot(x_mod,polymod(x_mod),'r', linestyle='dashed',alpha=0.5)
                else:
                    plt.plot(x, intercept + (slope*np.array(x)), 'r', linestyle='dashed',alpha=0.5)
                plt.title('r = '+str(np.round(r_value,2))+' p ='+str(np.round(p_value,2)))

              ax.spines['right'].set_visible(False)
              ax.spines['top'].set_visible(False)
              plt.tight_layout()
              return x,y



          ##############################
          #ERP
          ##############################
          def calc_erp(self,ev_type = None, ev_item = None,time_bin_size_ms = 100,apply_gauss_smoothing = True,gauss_sd_scaling = 0.1):
              # calculates ERP. (automatically calculates raw and z-scored). Calculates an absolute(z-scored ERP binned into time bins) that can be used to aggregate data across trials and subjects
              # Inputs:
              # ev_type...'WAIT','INSTRUCT', or 'MOVE'; if None, includes all events
              # ev_item...'Left Hand', 'Mouth and Tongue', 'Right Hand'; if None, includes all events
              #time_bin_size_ms = 100 ... this sets time bin length for time window binning 
              # apply_gauss_smoothing = False ..... if true, applies gaussian smoothing to abs_z_erp_ms ONLY (intended for visualization of aggregate data)
              # gauss_sd_scaling = 0.1.. (samplerate * 0.1) sets the width of gaussian to smooth with.



              # get eeg array
              eeg = self.elec_dict['subjEEG']

              # z-score over entire session (all ev types)
              eeg_mean = np.nanmean(eeg)
              eeg_sd = np.nanstd(eeg)

              # filter by event type
              if ev_type == None:
                         type_idx = np.ones((len(self.elec_dict['ev_type']))).astype('bool')
              else:
                         type_idx = self.elec_dict['ev_type'] == ev_type

              # filter by event item
              if ev_item == None:
                         item_idx = np.ones((len(self.elec_dict['ev_item']))).astype('bool')
              else:
                         item_idx = self.elec_dict['ev_item'] == ev_item

              #calculate erp
              erp = eeg[type_idx&item_idx,:]
              erp_mean = np.nanmean(erp,axis=0)
              erp_sem = stats.sem(erp,axis=0)

              #drop buffer and update self.elec_dict
              self.elec_dict['erp'] = self.dropBuffer(erp)
              self.elec_dict['erp_mean'] = self.dropBuffer(erp_mean)
              self.elec_dict['erp_sem'] = self.dropBuffer(erp_sem)

              # calculate xval vector (for plotting); this assumes buffer has been removed
              self.elec_dict['erp_xval'] = np.linspace(self.ms_to_samples(self.params['tmin_ms']), self.ms_to_samples(self.params['tmax_ms']),len(self.elec_dict['erp_mean']))

              # save params
              self.elec_dict['erp_params'] = {'ev_type':ev_type,'ev_item':ev_item,'time_bin_size_ms':time_bin_size_ms,'apply_gauss_smoothing':apply_gauss_smoothing,'gauss_sd_scaling':gauss_sd_scaling}

              # z scored eeg mean
              self.elec_dict['z_erp'] = (self.elec_dict['erp'] - eeg_mean)/(eeg_sd)
              self.elec_dict['abs_z_erp_mean'] = np.absolute(np.nanmean(self.elec_dict['z_erp'],axis=0))
              self.elec_dict['z_erp_sem'] = stats.sem(self.elec_dict['z_erp'],axis=0)

              # resample in ms (so you can average across patients)
              self.elec_dict['abs_z_erp_mean_ms'] = self.data2timeBins(self.elec_dict['abs_z_erp_mean'],time_bin_size_ms=time_bin_size_ms)
              self.elec_dict['z_erp_sem_ms'] = self.data2timeBins(self.elec_dict['z_erp_sem'],time_bin_size_ms=time_bin_size_ms)              
              self.elec_dict['erp_xval_ms'] = self.time_bin_dict['lbl']


              if apply_gauss_smoothing == True:
                         self.elec_dict['abs_z_erp_mean_ms'] =ndimage.filters.gaussian_filter1d(self.elec_dict['abs_z_erp_mean_ms'],sigma =(gauss_sd_scaling*self.elec_dict['sr']))


          def plot_erp(self,ax = None,figsize=None,fsize_tick=12,alpha=0.5,xRange_ms = (-500,1000),apply_time_bins = False,yL=None):
              # assumes that data is in samples
              if ax == None:
                         f = plt.figure(figsize=figsize)              
                         ax = plt.subplot(111)

              # get xval, mean and sem
              if apply_time_bins == False:
              	xval = self.elec_dict['erp_xval']
              	erp_mean = self.elec_dict['erp_mean']
              	erp_sem = self.elec_dict['erp_sem']
              else:
              	# use binned data
              	xval = self.elec_dict['erp_xval_ms']
              	erp_mean = self.elec_dict['erp_mean_ms']
              	erp_sem = self.elec_dict['erp_sem_ms']


              # plot it
              plt.plot(xval,erp_mean,alpha=alpha,label = self.elec_dict['erp_params']['ev_type']+'/'+str(self.elec_dict['erp_params']['ev_item']))
              plt.fill_between(xval,erp_mean+erp_sem,erp_mean-erp_sem,alpha=0.5)

              # find xRange indices
              xRange_samp_start = self.ms_to_samples(xRange_ms[0])
              xRange_samp_end = self.ms_to_samples(xRange_ms[1])

              # set x ticks
              if apply_time_bins == False:
              	# use sample data
              	xt = np.array([xRange_samp_start,0,0.5*self.elec_dict['sr'],xRange_samp_end])
              	ax.set_xticks(xt)
              	ax.set_xticklabels(np.round(1000*(xt/self.elec_dict['sr'])).astype('int'),fontsize=fsize_tick)
              	#set x lims
              	ax.set_xlim(xRange_samp_start,xRange_samp_end)
              else:
              	# use time bins data
                         xt = np.array([xRange_ms[0],0,xval[np.argmin(np.abs(xval-500))],xRange_ms[1]])
                         ax.set_xticks(xt)
                         ax.set_xticklabels(xt,fontsize=fsize_tick)
                         ax.set_xlim(xRange_ms[0],xRange_ms[1])

              if yL!=None:
                         ax.set_ylim(yL[0],yL[1])


          def plot_erp_byType(self,ev_item=None,apply_time_bins=False,apply_gauss_smoothing = True,gauss_sd_scaling = 0.1,yL=None):
              # if ev_item = None, it will combine data from all items. Otherwise, it will only plot data from the specified item (left/right/mouth)
              f = plt.figure()
              ax = plt.gca()
              evType_list = np.unique(self.elec_dict['ev_type'])

              for t in evType_list:
                         self.calc_erp(ev_type=t,ev_item = ev_item,apply_gauss_smoothing = apply_gauss_smoothing,gauss_sd_scaling = gauss_sd_scaling)
                         self.plot_erp(ax=ax,apply_time_bins=apply_time_bins,yL=yL)

              plt.legend()
              plt.vlines(x=0,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],alpha = 0.3,color='r')
              plt.xlabel('Time from stimulus (ms)')
              plt.ylabel('Voltage (uV)')
              

          def plot_erp_byItem(self,ev_type = 'MOVE',apply_time_bins=False,apply_gauss_smoothing = True,gauss_sd_scaling = 0.1,yL=None):

              f = plt.figure()
              ax = plt.gca()

              evItem_list = np.unique(self.elec_dict['ev_item'])

              for i in evItem_list:
                         self.calc_erp(ev_type=ev_type,ev_item=i,apply_gauss_smoothing = apply_gauss_smoothing,gauss_sd_scaling = gauss_sd_scaling)
                         self.plot_erp(ax=ax,apply_time_bins=apply_time_bins,yL=yL)                                  
              plt.legend()
              plt.vlines(x=0,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],alpha = 0.3,color='r')
              plt.xlabel('Time from stimulus (ms)')
              plt.ylabel('Voltage (uV)')

          ##############################
          #LFO
          ##############################
          # find LFO
          def calc_lfo(self,ev_type=None,ev_item = None,offset_ms = 0, duration_ms = 5000,fRange = None,aperiodic_mode='knee'):
              # this function applies FOOF on continous EEG data
              # it updates params dict with frange_LFO which is the frequency range of
              # the dominant low freq oscillation (LFO)
              # input:
              # eeg_vec_w_buff 
              # parse fRange
              if fRange == None:
                         fRange = self.params['frange_LFA']

              #get eeg_vec
              eeg_vec_w_buff = self.catEEG(ev_type=ev_type,ev_item=None,offset_ms = offset_ms, duration_ms = duration_ms)

              # calculate PSD
              psd,freqs = mne.time_frequency.psd_array_welch(eeg_vec_w_buff,sfreq = self.elec_dict['sr'], fmin=fRange[0],fmax=fRange[1])

              # initialize fooof model
              fm = ff.FOOOF(aperiodic_mode=aperiodic_mode) # using knee parameter leads to a better fit over long freq ranges

              # initialize dict to save results
              self.lfo_dict = self.initResultsDict()

              # Fit FOOOF model
              fm.fit(freqs,psd,fRange)

              # save data
              # full model 
              self.lfo_dict['fm'] = fm

              #collect model fit results in a dictionary (aperiodic params, peak params, r sq, fit error, gaussian params)
              # aperiodic params are a tuple, which represent (offset, knee, exponent)
              # peak params are tuple that represent (center freq, power, and band width for each peak)
              
              # freqs,psd
              self.lfo_dict['freqs'] = freqs
              self.lfo_dict['psd'] = psd  

              self.lfo_dict['aperiodic'],self.lfo_dict['peaks'],self.lfo_dict['rsq'],self.lfo_dict['fit_error'],self.lfo_dict['gauss_params'] = fm.get_results()

              # num peaks
              self.lfo_dict['num_peaks'] = np.shape(self.lfo_dict['peaks'])[0]

              # if there are peaks in the PSD
              if len(self.lfo_dict['peaks'])==0:
                         self.lfo_dict['LFO'] = (np.nan,np.nan,np.nan)
                         self.lfo_dict['LFO_fRange'] = (np.nan,np.nan)
                         self.updateParamsDict('frange_LFO',(np.nan,np.nan))
              else:
                         # find LFO (lowest freq peak to estimate ~theta-alpha range effects)

                         lfo_idx = np.argmin(self.lfo_dict['peaks'][:,0])  
                         self.lfo_dict['LFO'] = self.lfo_dict['peaks'][lfo_idx]

                         # get custom low frequency freq range
                         LFO_fRange= ((self.lfo_dict['LFO'][0]-(self.lfo_dict['LFO'][2]/2)),(self.lfo_dict['LFO'][0]+(self.lfo_dict['LFO'][2]/2)))

                         # make sure the floor is set to fRange[0]
                         if LFO_fRange[0]<2:
                                  LFO_fRange= (fRange[0],LFO_fRange[1])
                                  
                         # update params dict with LFO frange (data-driven freq range for this electrode). Can be usef for subsequent analyses
                         self.lfo_dict['LFO_fRange'] = LFO_fRange
                         self.updateParamsDict('frange_LFO',np.copy(self.lfo_dict['LFO_fRange']))
              
          
              # get THETA peak (within fRange_theta)
              #tuple representing center freq, power, band width
              self.lfo_dict['theta'] = ff.analysis.get_band_peak_fm(fm,
              self.params['frange_theta'], select_highest=True)


              # get alpha peak (within fRange_theta and alpha)
              #tuple representing center freq, power, band width
              self.lfo_dict['alpha'] = ff.analysis.get_band_peak_fm(fm,
              self.params['frange_alpha'], select_highest=True)          


              # get beta peak (within fRange_theta and alpha)
              #tuple representing center freq, power, band width
              self.lfo_dict['beta'] = ff.analysis.get_band_peak_fm(fm,
              self.params['frange_beta'], select_highest=True)              

              # get gamma peak (within fRange_theta and alpha)
              #tuple representing center freq, power, band width
              self.lfo_dict['gamma'] = ff.analysis.get_band_peak_fm(fm,
              self.params['frange_gamma'], select_highest=True)                        

              # store params
              self.lfo_dict['params'] = {'ev_type':ev_type,'offset_ms':offset_ms,'duration_ms':duration_ms,'fRange':fRange,'aperiodic_mode':aperiodic_mode}

          def plot_lfo(self, ax = None):
              # plots cached fooof data
              if ax == None:
                         f = plt.figure()
                         ax = plt.subplot(111)

              self.lfo_dict['fm'].print_results()
              self.lfo_dict['fm'].plot(plot_peaks ='line',ax=ax)

          
          ##############################
          #SYNCHRONY
          ##############################
          # Two options on how to compute synchrony (methods 2 and 3 are intended for use with phase difference distributions with other)

          # 1) method = 'xTrial'. Here it computes across trial phase consistency (on sample and time binned data). It generates a null distribution using a permutation perocedure on the entire phase matrix (or phase difference matrix). Randomizes the order of trials and circularly shifts each trial. 

          # 2) method = 'xTime'. Here it computes within trial phase consistency and aggregates r-square values across trials (on sample and time binned data). It generates a null distribution for phase reset by randomly sampling the entire phase matrix. 

          # Note: I considered computing 'xSess' measures that involve concatonating the entire session into a single time series, however, these are not appropriate because the value of the phase difference may change throughout the session resulting in a uniform distribution over a long enough time scale

          def calc_phaseSync_xTrial(self, other = None, fRange=(3,8),ev_type = None, ev_item = None,time_bin_size_ms = 100, method = 'filt_hilb',wave_numfreqs = None, wave_numcycles = None,n_iters=10000):
              # this computes across trial phase consistency (on sample and time binned data).  if other is another electrode object, it will calculate across trial consistency of phase differences [] Need to implement NULL - For null generates a null distribution using a permutation perocedure on the entire phase matrix (or phase difference matrix). Randomizes the order of trials and circularly shifts each trial. 

              # load phase
              phsmat_trFilt = self.get_phase_filteredByTrial(other = other, fRange = fRange,ev_type = ev_type,ev_item = ev_item, method = method,wave_numfreqs = wave_numfreqs, wave_numcycles = wave_numcycles)

              # bin into time bins (will work on this in parallel)
              phsmat_trFilt_ms = self.data2timeBins(phsmat_trFilt,time_bin_size_ms = time_bin_size_ms, isCircData = True)


              # initialize phs dict for results
              self.phsSyncDict_xTrial = self.initResultsDict(other = other,incl_LFO=True)


              # populate it
              #x-val data
              self.phsSyncDict_xTrial['xval'] = np.linspace(self.ms_to_samples(self.params['tmin_ms']),self.ms_to_samples(self.params['tmax_ms']), np.shape(phsmat_trFilt)[1])

              self.phsSyncDict_xTrial['xval_ms'] = self.time_bin_dict['lbl']



              self.phsSyncDict_xTrial['params'] ={'ev_type':ev_type,'ev_item':ev_item,'fRange':fRange,'time_bin_size_ms':time_bin_size_ms}


              # RVL - resultant vector length
              self.phsSyncDict_xTrial['RVL'] = circ.resultant_vector_length(phsmat_trFilt,axis=0)
              self.phsSyncDict_xTrial['RVL_bin_ms'] = circ.resultant_vector_length(phsmat_trFilt_ms,axis=0)  


              # ral z test (assumes unimodal data). Not using omnibus test here because there may be too few trials
              self.phsSyncDict_xTrial['ral_p'],self.phsSyncDict_xTrial['ral_z'] = circ.rayleigh(phsmat_trFilt,axis=0)
              self.phsSyncDict_xTrial['ral_p_ms'],self.phsSyncDict_xTrial['ral_z_ms'] = circ.rayleigh(phsmat_trFilt_ms,axis=0)

              # preferred phase
              self.phsSyncDict_xTrial['prefPhase']  = circ.mean(phsmat_trFilt,axis=0)
              self.phsSyncDict_xTrial['prefPhase_ms']  = circ.mean(phsmat_trFilt_ms,axis=0)

              #implement Null dist here
              nullphs = self.calc_null_phase(n_iters = n_iters)

              # bin into time bins
              nullphs_ms = self.data2timeBins(nullphs,time_bin_size_ms = time_bin_size_ms, isCircData = True)

              # Null stats
              self.phsSyncDict_xTrial['RVL_null'] = circ.resultant_vector_length(nullphs,axis=0)
              self.phsSyncDict_xTrial['RVL_null_bin_ms'] = circ.resultant_vector_length(nullphs_ms,axis=0)  


          def calc_phaseSync_xTime(self, other = None, fRange=(3,8),ev_type = None, ev_item = None, method = 'filt_hilb',wave_numfreqs = None, wave_numcycles = None,offset_ms_list=[0], duration_ms=1000,n_iters=10000):
                #Here we computes within trial phase consistency. Several measures of phase consistency (RVL (assumes unimodal),and omnibus (less sensitive, but works with multimodal as well); it aggregates these measures across trials. It also computes xSess measures on time series data  of across session phase consistency by concatenating phase differences across trials into a single time series
                # inputs:
                # other ... another electrode object. Cannot be None because this mehtod in not appropriate for single trial analysis
                # fRange ... freq range to focus on (customize this based on the LFO at these electrode pairs)
                # ev_item ... option to filter by items
                # ev_type .. option to filter by ev type 
                # method ... 'filt_hilb' vs 'wave'
                # wave_numfreqs and wave_numcycles are parameters for wavelet analyses (if method = wave)
                # offset_ms_list = [0] ... list of start times relative to event onset (can range from params['tmin'] to params['tmax'])

                # error check:
                if other==None:
                     raise NameError('OTHER must be another Electrode object. This method is not appropriate for single electrode analyses')

                # load phase difference data 
                phsmat_trFilt = self.get_phase_filteredByTrial(other = other, fRange = fRange,ev_type = ev_type,ev_item = ev_item, method = method,wave_numfreqs = wave_numfreqs, wave_numcycles = wave_numcycles)

                # GET NULL data

                #load null phase mat (based on params; dont need to filter by trial). Clip based on duration_ms (use this distribution to perform non-parametric stats below)
                nullphs = self.calc_null_phase(n_iters = n_iters)

                # clip null phase data based on duration_ms
                nullphs_clip,nullxval = self.clip_data(datamat_nobuff=nullphs, offset_ms=0, duration_ms=duration_ms)

                # calc null distribution of statistics
                null_RVL = circ.resultant_vector_length(nullphs_clip,axis=1)
                null_RVL_z = (null_RVL-np.nanmean(null_RVL))/np.nanstd(null_RVL)
                null_omn_p, null_omn_m = circ.omnibus(nullphs_clip,axis=1)


                # initialize phs dict for results
                self.phsSyncDict_xTime = self.initResultsDict(other = other,incl_LFO=True)

                # store params
                self.phsSyncDict_xTime['params'] ={'ev_type':ev_type,'ev_item':ev_item,'fRange':fRange,'offset_ms':offset_ms_list,'duration_ms':duration_ms}

                # initialize fields to populate during the loop
                # lists  - dont need to concatonate at Subject and Group level
                self.phsSyncDict_xTime['phsmat_clip_list'] = []
                self.phsSyncDict_xTime['xvalms_clip_list'] = []
                
                # arrays - will be helpful for concatonattion later at Subj and Group levels
                self.phsSyncDict_xTime['xval_ms'] = np.empty(0) # this can be used to plot
                self.phsSyncDict_xTime['RVL']  = np.empty(0)
                self.phsSyncDict_xTime['RVL_null']  = null_RVL 
                self.phsSyncDict_xTime['RVL_mean'] = np.empty(0)
                self.phsSyncDict_xTime['RVL_mean_null'] = np.empty(0)
                self.phsSyncDict_xTime['RVL_mean_zVsNull'] = np.empty(0)
                self.phsSyncDict_xTime['RVL_mean_zVsNull_null'] = np.empty(0)

                self.phsSyncDict_xTime['RVL_tVsNull_t'] = np.empty(0)
                self.phsSyncDict_xTime['RVL_tVsNull_p'] = np.empty(0)

                # omnibus test
                self.phsSyncDict_xTime['omn_p'] = np.empty(0)
                self.phsSyncDict_xTime['omn_p_vsNull'] = np.empty(0)  
                self.phsSyncDict_xTime['omn_z_vsNull'] = np.empty(0)
                self.phsSyncDict_xTime['omn_tVsNull_t'] = np.empty(0)
                self.phsSyncDict_xTime['omn_tVsNull_p'] = np.empty(0)

                # pref phase
                self.phsSyncDict_xTime['prefPhase']  =  np.empty(0)

                # Circ ANOVA - ev_item [] can compute this on z-vals for now (to assess whether there is more consistency for particular movements or not)
                self.phsSyncDict_xTime['circanov_evitem_p']  =  np.empty(0)

                # loop through time windows and perform computations
                for t in offset_ms_list:
                     phsmat_clip,xval_ms_clip = self.clip_data(datamat_nobuff=phsmat_trFilt, offset_ms=t, duration_ms=duration_ms)

                     # calculate measures
                     RVL = circ.resultant_vector_length(phsmat_clip,axis=1)
                     omn_p, omn_m = circ.omnibus(phsmat_clip,axis=1)
                     pref_phase = circ.mean(phsmat_clip,axis=1)

                     # compare these with null
                     RVL_zByNull = (RVL-np.nanmean(null_RVL))/np.nanstd(null_RVL)


                     # calculate non-parametric omnibus p-values by comparing to null distribution
                     omn_p_vsNull = np.empty(0)
                     for i in np.arange(0,len(omn_p)):
                          omn_p_vsNull=np.append(omn_p_vsNull,np.count_nonzero(null_omn_p <= omn_p[i])/len(null_omn_p))

                     # transform to z-score: first get probability of false positive as (1-non-param p). Then use stats.norm.ppf to z-score. Such that non-param p = 0.05 -> prob of 0.95 -> z = 1.96
                     omn_z_vsNull = stats.norm.ppf(1-omn_p_vsNull)


                     # populate results dict

                     # lists
                     self.phsSyncDict_xTime['phsmat_clip_list'].append(phsmat_clip)
                     self.phsSyncDict_xTime['xvalms_clip_list'].append(xval_ms_clip)

                     # 1d arrays
                     self.phsSyncDict_xTime['xval_ms'] = np.append(self.phsSyncDict_xTime['xval_ms'],np.round(np.nanmean(xval_ms_clip))) # for plotting
                     # RVL mean
                     self.phsSyncDict_xTime['RVL_mean'] = np.append(self.phsSyncDict_xTime['RVL_mean'],np.nanmean(RVL)) 
                     self.phsSyncDict_xTime['RVL_mean_null'] = np.append(self.phsSyncDict_xTime['RVL_mean_null'],np.nanmean(null_RVL)) 

                     self.phsSyncDict_xTime['RVL_mean_zVsNull'] = np.append(self.phsSyncDict_xTime['RVL_mean_zVsNull'],np.nanmean(RVL_zByNull)) 
                     self.phsSyncDict_xTime['RVL_mean_zVsNull_null'] = np.append(self.phsSyncDict_xTime['RVL_mean_zVsNull_null'],np.nanmean(null_RVL_z))
                    

                     # RVL_tVsNull .. overall measure of phase synch for this clip across all trials
                     RVL_tVsNull_t,RVL_tVsNull_p= stats.ttest_1samp(RVL_zByNull,popmean=0)
                     self.phsSyncDict_xTime['RVL_tVsNull_t'] = np.append(self.phsSyncDict_xTime['RVL_tVsNull_t'],RVL_tVsNull_t)
                     self.phsSyncDict_xTime['RVL_tVsNull_p'] = np.append(self.phsSyncDict_xTime['RVL_tVsNull_p'],RVL_tVsNull_p)

                     # omn_tVsNull... overall measure of phase synch for this clup across trials (does not assume unimodal data)
                     omn_tVsNull_t,omn_tVsNull_p= stats.ttest_1samp(omn_z_vsNull,popmean=0)
                     self.phsSyncDict_xTime['omn_tVsNull_t'] = np.append(self.phsSyncDict_xTime['omn_tVsNull_t'],omn_tVsNull_t)
                     self.phsSyncDict_xTime['omn_tVsNull_p'] = np.append(self.phsSyncDict_xTime['omn_tVsNull_p'],omn_tVsNull_p)

                     # 2d arrays
                     # RVL by trial
                     if t == offset_ms_list[0]:
                          # first iter
                          #RVL (raw,null and z-scored vs. null)
                          self.phsSyncDict_xTime['RVL'] = RVL[:,np.newaxis]
                          self.phsSyncDict_xTime['RVL_zByNull'] = RVL_zByNull[:,np.newaxis]

                          # omnibus
                          self.phsSyncDict_xTime['omn_p'] = omn_p[:,np.newaxis]
                          self.phsSyncDict_xTime['omn_p_vsNull'] = omn_p_vsNull[:,np.newaxis]
                          self.phsSyncDict_xTime['omn_z_vsNull'] = omn_z_vsNull[:,np.newaxis]

                          # pref phase
                          self.phsSyncDict_xTime['prefPhase']  = pref_phase[:,np.newaxis]
                     else:
                          # RVL (raw and z by null)
                          self.phsSyncDict_xTime['RVL'] = np.append(self.phsSyncDict_xTime['RVL'],RVL[:,np.newaxis],axis =1)
                          self.phsSyncDict_xTime['RVL_zByNull'] = np.append(self.phsSyncDict_xTime['RVL_zByNull'],RVL_zByNull[:,np.newaxis],axis =1)

                          # omnibus (raw and vs null)
                          self.phsSyncDict_xTime['omn_p'] = np.append(self.phsSyncDict_xTime['omn_p'],omn_p[:,np.newaxis],axis =1)
                          self.phsSyncDict_xTime['omn_p_vsNull'] = np.append(self.phsSyncDict_xTime['omn_p_vsNull'],omn_p_vsNull[:,np.newaxis],axis =1)
                          self.phsSyncDict_xTime['omn_z_vsNull'] = np.append(self.phsSyncDict_xTime['omn_z_vsNull'],omn_z_vsNull[:,np.newaxis],axis =1)

                          # pref phase
                          self.phsSyncDict_xTime['prefPhase']  = np.append(self.phsSyncDict_xTime['prefPhase'],pref_phase[:,np.newaxis],axis =1)

                     # first pass analysis: One way anova to assess item coding
                     if (ev_item == None):
                          if ev_type != None:
                                ev_item_list = self.phsmat_dict['ev_item'][self.phsmat_dict['ev_type'] == ev_type]
                          else:
                                ev_item_list = self.phsmat_dict['ev_item']

                          f,p = stats.f_oneway(omn_z_vsNull[ev_item_list=='Right Hand'] ,omn_z_vsNull[ev_item_list=='Left Hand'] ,omn_z_vsNull[ev_item_list=='Mouth and Tongue'])

                          self.phsSyncDict_xTime['circanov_evitem_p'] = np.append(self.phsSyncDict_xTime['circanov_evitem_p'],p)
                     else:
                          self.phsSyncDict_xTime['circanov_evitem_p'] = np.append(self.phsSyncDict_xTime['circanov_evitem_p'],np.nan)



          #https://stackoverflow.com/questions/22562364/circular-histogram-for-python
          def plot_phaseSync_rose(self,ax, angles, bins=16, density=None, offset=0, lab_unit="degrees",start_zero=False,ecol='0.5', **param_dict):
            """
            Plot polar histogram of angles on ax. ax must have been created using
            subplot_kw=dict(projection='polar'). Angles are expected in radians.
            """
            # Wrap angles to [-pi, pi)
            angles = (angles + np.pi) % (2*np.pi) - np.pi

            # Set bins symetrically around zero
            if start_zero:
                 # To have a bin edge at zero use an even number of bins
                 if bins % 2:
                      bins += 1
                 bins = np.linspace(-np.pi, np.pi, num=bins+1)

            # Bin data and record counts
            count, bin = np.histogram(angles, bins=bins)

            # Compute width of each bin
            widths = np.diff(bin)

            # By default plot density (frequency potentially misleading)
            if density is None or density is True:
                 # Area to assign each bin
                 area = count / angles.size
                 # Calculate corresponding bin radius
                 radius = (area / np.pi)**.5
            else:
                 radius = count

            # Plot data on ax
            ax.bar(bin[:-1], radius, zorder=1, align='edge', width=widths,
                     edgecolor=ecol, fill=False, linewidth=1)

            # Set the direction of the zero angle
            ax.set_theta_offset(offset)

            # Remove ylabels, they are mostly obstructive and not informative
            ax.set_yticks([])

            if lab_unit == "radians":
                 label = ['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$',
                              r'$\pi$', r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$']
                 ax.set_xticklabels(label)

          def plot_phaseSync_xTime_rvlByTime(self,offset_of_interest_ms = 100,ax = None,yL = None):
            # This function plots RVL over time for phase synchrony results that are cached in phsSyncDict_xTime. It marks the null distribution and a time window of interest

          
            # plot RVL_mean over time indicating RVL_mean_null
            if ax == None:
                f = plt.figure()
                ax = plt.subplot(111)
            ax.plot(self.phsSyncDict_xTime['xval_ms'],self.phsSyncDict_xTime['RVL_tVsNull_t'])    
            #ax.plot(self.phsSyncDict_xTime['xval_ms'],self.phsSyncDict_xTime['RVL_mean'])
            #ax.plot(self.phsSyncDict_xTime['xval_ms'],self.phsSyncDict_xTime['RVL_mean_null'],color = '0.5',linestyle = 'dashed')
            # option to hard set the y lim
            if yL != None:
                  ax.set_ylim(yL)

            # find the time window of interest per input param
            t_idx = np.argmin(np.abs(self.phsSyncDict_xTime['params']['offset_ms'] - offset_of_interest_ms))
            ax.vlines(self.phsSyncDict_xTime['xval_ms'][t_idx],ax.get_ylim()[0],ax.get_ylim()[1],linestyle='solid',color = 'r',alpha=0.5)

            # set labels
            ax.set_xlabel('Time from cue (ms), duration of each time window = '+str(self.phsSyncDict_xTime['params']['duration_ms']))
            ax.set_ylabel('Synchrony t-statistic')
            ax.axis('tight')


          def plot_phaseSync_xTime_rvlVsNull(self,offset_of_interest_ms = 100,ax=None):
            # This function plots RVL for each trial vs. the null RVL distribution for a given time window
          
            if ax == None:
                plt.figure()
                ax = plt.subplot(111)  

            # find time window of interest
            t_idx = np.argmin(np.abs(self.phsSyncDict_xTime['params']['offset_ms'] - offset_of_interest_ms))


            rvl_true = self.phsSyncDict_xTime['RVL'][:,t_idx]
            rvl_null = self.phsSyncDict_xTime['RVL_null']
            ax.hist(rvl_null,bins=100,edgecolor=None,alpha=0.5,color='0.6')
            ax.vlines(rvl_true,ax.get_ylim()[0],ax.get_ylim()[1],linestyle='solid',color = 'r',alpha=0.2)
            ax.vlines(np.nanmean(rvl_true),ax.get_ylim()[0],ax.get_ylim()[1],linestyle='solid',color = 'b',alpha=1)
            ax.set_title('t vs. null = '+str(np.round(self.phsSyncDict_xTime['RVL_tVsNull_t'][t_idx],2))+', p = '+str(np.round(self.phsSyncDict_xTime['RVL_tVsNull_p'][t_idx],4)))
            ax.set_xlabel('Resultant vector length')
            ax.set_ylabel('count')




          def plot_phaseSync_xTime_eegClips(self,offset_of_interest_ms = 100,n_clips = 3,n_phs_bins=30,savefig=None,figdir = None, fname=None):
            # This function plots eeg clips of raw and filtered eeg for the highest synchrony trials. Creates a new figure for each clip


            # find time window of interest
            t_idx = np.argmin(np.abs(self.phsSyncDict_xTime['params']['offset_ms'] - offset_of_interest_ms))
            rvl_true = self.phsSyncDict_xTime['RVL'][:,t_idx]

            #sort trials from high to lowsynchrony
            sync_idx = np.flip(np.argsort(rvl_true)) 

            # get eeg clips  
            offset_ms = offset_of_interest_ms
            duration_ms = self.phsSyncDict_xTime['params']['duration_ms']
            fRange = self.phsSyncDict_xTime['params']['fRange']
            eeg_clip_a,xval_ms = self.clip_data(self.dropBuffer(self.phsmat_dict['eeg_self']),offset_ms=offset_ms,duration_ms=duration_ms)
            eeg_clip_b,xval_ms = self.clip_data(self.dropBuffer(self.phsmat_dict['eeg_other']),offset_ms=offset_ms,duration_ms=duration_ms)

            # filter eeg clips 
            eeg_clip_a_filt = mne.filter.filter_data(data = eeg_clip_a.astype('float'),
                                                                sfreq = self.elec_dict['sr'],
                                                                l_freq = fRange[0],
                                                                h_freq = fRange[1],
                                                                method='fir',verbose=False)
            eeg_clip_b_filt = mne.filter.filter_data(data = eeg_clip_b.astype('float'),
                                                                sfreq = self.elec_dict['sr'],
                                                                l_freq = fRange[0],
                                                                h_freq = fRange[1],
                                                                method='fir',verbose=False)

            for i in np.arange(0,n_clips):
                phsmat_high_sync = self.phsSyncDict_xTime['phsmat_clip_list'][t_idx][sync_idx[i]]
                f = plt.figure()
                gs = f.add_gridspec(nrows=2,ncols=3)

                # plot polar histogram of phase differences for this clip
                ax = plt.subplot(gs[0,0],projection='polar')
                #ax=plt.hist(phsmat_high_sync,bins=50)
                self.plot_phaseSync_rose(ax=ax,angles= phsmat_high_sync, bins=n_phs_bins, density=None, offset=0, lab_unit="radians",ecol='C0')
                

                ax2 = plt.subplot(gs[1,0],projection='polar')
                ax2.plot([0,circ.mean(phsmat_high_sync)],[0,rvl_true[sync_idx[i]]],linewidth=3,color = 'C1',alpha = 0.5)
                ax2.set_rlim([0,1])
                ax2.set_rticks([])
                #ax2.set_xticks([np.deg2rad(0),np.deg2rad(90),np.deg2rad(180),np.deg2rad(270)])
                #ax2.set_xticklabels(labels=[0,90,180,270])
                ax2.set_xticks([])
                ax2.text(np.deg2rad(200),.75,'RVL = '+str(np.round(rvl_true[sync_idx[i]],2)),fontsize=12)
                #ax2.set_xticklabels(labels=['$0$', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$',r'$\pi$', r'$5\pi/4$', r'$3\pi/2$', r'$7\pi/4$'])
                ax2.grid(False)

                


                # plot raw eeg
                ax = plt.subplot(gs[0,1:])
                plt.plot(xval_ms,eeg_clip_a[sync_idx[i],:],color='red',linestyle='solid',alpha=0.5)
                plt.plot(xval_ms,eeg_clip_b[sync_idx[i],:],color='blue',linestyle='solid',alpha=0.5)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.set_title('Raw EEG',fontsize=20)


                ax = plt.subplot(gs[1,1:])
                plt.plot(xval_ms,eeg_clip_a_filt[sync_idx[i],:],color='red',linestyle='dashed',alpha=0.5)
                plt.plot(xval_ms,eeg_clip_b_filt[sync_idx[i],:],color='blue',linestyle='dashed',alpha=0.5)
                ax.spines['left'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.set_title('Filtered EEG '+str(fRange[0])+'-'+str(fRange[1])+' Hz',fontsize=20)
                ax.set_yticks([])
                ax.set_xlabel('Time from cue (ms)')


                # save fig


                if savefig ==True:
                    ax = plt.gca()
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)
                    plt.tight_layout()
                    plt.savefig(fname = figdir+fname+'_eegClip'+str(i)+'.pdf')


          def plot_phaseSync_xTime_summary(self,offset_of_interest_ms = 100,n_clips = 3,n_phs_bins=30,ax_rvlByTime =None,ax_rvlVsNull=None,savefig = False, figdir=None,fname=None,fsize_ticks = 14,fsize_lbl = 20,yL = None):

            # plot rvl by time
            self.plot_phaseSync_xTime_rvlByTime(offset_of_interest_ms = offset_of_interest_ms,ax=ax_rvlByTime,yL=yL)

            # save fig 
            if savefig ==True:
                ax = plt.gca()
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.set_xlabel('Time from cue (ms)',fontsize = fsize_lbl)
                ax.set_ylabel('Synchrony t-statistic',fontsize = fsize_lbl)
                ax.set_xticklabels(ax.get_xticks().astype('int'),fontsize=fsize_ticks)
                ax.set_yticklabels(np.round(ax.get_yticks(),2),fontsize=fsize_ticks)
                ax.set_title('')
                plt.tight_layout()
                plt.savefig(fname = figdir+fname+'_rvlByTime.pdf')

            # plot rvl vs. null
            self.plot_phaseSync_xTime_rvlVsNull(offset_of_interest_ms = offset_of_interest_ms,ax=ax_rvlVsNull)
            # save fig 
            if savefig ==True:
                ax = plt.gca()
                ax.spines['right'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.set_xlabel('Null resultant vector length',fontsize = fsize_lbl)
                ax.set_ylabel('Count',fontsize = fsize_lbl)
                ax.set_xticklabels(np.round(ax.get_xticks(),2),fontsize=fsize_ticks)
                ax.set_yticklabels(ax.get_yticks().astype('int'),fontsize=fsize_ticks)
                print(ax.get_title())
                ax.set_title('')
                plt.tight_layout()
                plt.savefig(fname = figdir+fname+'_rvlVsNull.pdf')


            #plot eeg clips
            self.plot_phaseSync_xTime_eegClips(offset_of_interest_ms = offset_of_interest_ms,n_clips = n_clips,n_phs_bins=n_phs_bins,savefig = savefig, figdir=figdir,fname=fname)







##############################################################################
#SUBJECT CLASS
###############################################################################

class Subject(Electrode):
          # Performs group-level analysis within a single subject (e.g., pairwise synchrony). Inherits from Electrode class and can instantiate multiple Electrode objects to perform self-other analyses
          # ti will collect electrode-level results in "subjDf_" fields and collapsed results for the subject in "subjDataCollapsed_" fields (as a dictionary)


          # constructor
          def __init__(self,subj,clusId_list = None,roi_list = None, paramsDict = None):
              # clusId_list... list of clusId values to include (e.g. [1,2,3] includes electrodes that were assigned to the first three clusters)
              # paramsDict .. this is optional and will overwrite default params in electrode class

              self.subj = subj
              self.clusId_list = clusId_list
              self.roi_list = roi_list
              self.paramsDict = paramsDict

              # get bool idx for subj electrodes
              idx_thisSubj = groupDf.eval('subj==@subj').to_numpy()

              # get bool idx for cluster electrodes
              if clusId_list==None:
                         # if none is provided, we will include all clusters
                         idx_thisClus = np.ones(len(idx_thisSubj)).astype('bool')
              else:
                         idx_thisClus = groupDf.eval('clusId in @clusId_list').to_numpy()

              # get bool idx for ROI_list
              if roi_list==None:
                         # if none is provided, we will include all clusters
                         idx_thisRegion = np.ones(len(idx_thisSubj)).astype('bool')
              else:
                         idx_thisRegion = groupDf.eval('ROI in @roi_list').to_numpy()

              # get idx for electrodes that belong to this subject and cluster list
              self.e_idx_thisSubjClus = (idx_thisSubj&idx_thisClus&idx_thisRegion).nonzero()[0]
              self.uElbls_thisSubjClus = groupDf['uElbl'].iloc[self.e_idx_thisSubjClus].to_numpy()
              if len(self.e_idx_thisSubjClus)==0:
                         print('No Matching electrodes')
                         self.isempty = True
                         return
              else:
                         self.isempty = False


              #initialize Electrode object with first electrode so we can access methods
              Electrode.__init__(self,uElbl = groupDf.iloc[self.e_idx_thisSubjClus[0]]['uElbl'],paramsDict=paramsDict)          
              

          ########################## 
          #GENERAL FUNCTIONS
          ########################## 

          #pairs of electrodes within group (self) and between groups (self,other)
          def make_pairs(self):
                # generates a pair of uElbls
                # create a dict to store data
                subjPair_dict = {}
                subjPair_dict['pair_uElbl'] = []
                subjPair_dict['pair_uElbl_self'] = []
                subjPair_dict['pair_uElbl_other'] = []

                #create a list of all possible pairs for this subject (based on clustId and ROI specifications)
                for i in np.arange(0,len(self.uElbls_thisSubjClus)):
                     for j in np.arange(i,len(self.uElbls_thisSubjClus)):
                          # skip if we find a pair of identical electrodes
                          if i == j:
                                continue
                          subjPair_dict['pair_uElbl'].append(self.uElbls_thisSubjClus[i]+'--'+self.uElbls_thisSubjClus[j])
                          subjPair_dict['pair_uElbl_self'].append(self.uElbls_thisSubjClus[i])
                          subjPair_dict['pair_uElbl_other'].append(self.uElbls_thisSubjClus[j])

                self.subjPair_dict = subjPair_dict

          ##############################
          #ERP
          ##############################

          def calc_erp(self,ev_type = None, ev_item = None,time_bin_size_ms = 100,apply_gauss_smoothing = True,gauss_sd_scaling = 0.1):
              # loops through electrodes and updates electrode dict. Subject class, overwrites Electrode.calc_erp

              # check if isempty
              if self.isempty == True:
                         print('No matching electrode')
                         return


              # initialize dict list to collect ERP data
              subjDictList_erp = []

              # key list
              keys_to_keep = ['clusId','sr','postInst_HFA', 'postMove_HFA', 'effectorSelective', 'instructionSelective', 'subj', 'eLbl', 'x', 'y', 'z', 'anat', 'ROI', 'uElbl', 'erp_mean', 'erp_xval','abs_z_erp_mean_ms','erp_xval_ms']

              #loop through electrode list
              for e in self.e_idx_thisSubjClus:
                         E = Electrode(uElbl = groupDf.iloc[e]['uElbl'],paramsDict=self.paramsDict)

                         # calcERP
                         E.calc_erp(ev_type = ev_type, ev_item = ev_item,time_bin_size_ms = time_bin_size_ms,apply_gauss_smoothing = apply_gauss_smoothing,gauss_sd_scaling = gauss_sd_scaling)

                         # get elect
                         elec_dict = E.elec_dict
                         all_keys = list(elec_dict.keys())

                         # remove keys 
                         for k in all_keys:
                                  if (k in keys_to_keep)==False:
                                            del elec_dict[k]

                         # update 
                         subjDictList_erp.append(E.elec_dict)


              # convert to dataframe
              self.subjDf_erp = pd.DataFrame(subjDictList_erp,index=self.uElbls_thisSubjClus)
              

              # collapse data within subject (as a dict)
              self.subjDataCollapsed_erp = {}
              self.subjDataCollapsed_erp['erp_mean'] = np.nanmean(np.stack(self.subjDf_erp['erp_mean'].to_numpy()),axis=0)
              self.subjDataCollapsed_erp['erp_sem'] = stats.sem(np.stack(self.subjDf_erp['erp_mean'].to_numpy()),axis=0)
              self.subjDataCollapsed_erp['erp_xval'] = self.subjDf_erp.iloc[0]['erp_xval']

              self.subjDataCollapsed_erp['abs_z_erp_mean_ms'] = np.nanmean(np.stack(self.subjDf_erp['abs_z_erp_mean_ms'].to_numpy()),axis=0)
              self.subjDataCollapsed_erp['z_erp_sem_ms'] = stats.sem(np.stack(self.subjDf_erp['abs_z_erp_mean_ms'].to_numpy()),axis=0)
              self.subjDataCollapsed_erp['erp_xval_ms'] = self.subjDf_erp.iloc[0]['erp_xval_ms']


              self.subjDataCollapsed_erp['erp_params'] = {'ev_type':ev_type,'ev_item':ev_item,'time_bin_size_ms':time_bin_size_ms}



          def plot_erp(self,ax = None,figsize=None,fsize_tick=12,alpha=0.5,xRange_ms = (-500,1000),apply_time_bins=False,yL=None):
          	# for class subject. overwrites Electrode.plot_erp
              # assumes that data is in samples
              if ax == None:
                         f = plt.figure(figsize=figsize)              
                         ax = plt.subplot(111)
              # get xval, mean and sem
              if apply_time_bins == False:
              	xval = self.subjDataCollapsed_erp['erp_xval']
              	erp_mean = self.subjDataCollapsed_erp['erp_mean']
              	erp_sem = self.subjDataCollapsed_erp['erp_sem']
              else:
              	# use binned data
              	xval = self.subjDataCollapsed_erp['erp_xval_ms']
              	erp_mean = self.subjDataCollapsed_erp['abs_z_erp_mean_ms']
              	erp_sem = self.subjDataCollapsed_erp['z_erp_sem_ms']


              # plot it
              plt.plot(xval,erp_mean,alpha=alpha,label = self.subjDataCollapsed_erp['erp_params']['ev_type']+'/'+str(self.subjDataCollapsed_erp['erp_params']['ev_item']))
              plt.fill_between(xval,erp_mean+erp_sem,erp_mean-erp_sem,alpha=0.5)

              # find xRange indices
              xRange_samp_start = self.ms_to_samples(xRange_ms[0])
              xRange_samp_end = self.ms_to_samples(xRange_ms[1])

              # set x ticks
              if apply_time_bins == False:
              	# use sample data
              	xt = np.array([xRange_samp_start,0,0.5*self.elec_dict['sr'],xRange_samp_end])
              	ax.set_xticks(xt)
              	ax.set_xticklabels(np.round(1000*(xt/self.elec_dict['sr'])).astype('int'),fontsize=fsize_tick)
              	#set x lims
              	ax.set_xlim(xRange_samp_start,xRange_samp_end)
              else:
              	# use time bins data
                         xt = np.array([xRange_ms[0],0,xval[np.argmin(np.abs(xval-500))],xRange_ms[1]])
                         ax.set_xticks(xt)
                         ax.set_xticklabels(xt,fontsize=fsize_tick)
                         ax.set_xlim(xRange_ms[0],xRange_ms[1])

              if yL!=None:
                         ax.set_ylim(yL[0],yL[1])

          ##############################
          #LFO
          ##############################

          def calc_lfo(self,ev_type = None, ev_item = None,offset_ms = 0, duration_ms = 1000):
              # loops through electrodes and updates electrode dict. Subject class, overwrites Electrode.calc_lfo

              # check if isempty
              if self.isempty == True:
                         print('No matching electrode')
                         return


              # initialize dict list to collect ERP data
              subjDictList_lfo = []

              #loop through electrode list
              for e in self.e_idx_thisSubjClus:
                         E = Electrode(uElbl = groupDf.iloc[e]['uElbl'],paramsDict=self.paramsDict)

                         # calcERP
                         E.calc_lfo(ev_type = ev_type, ev_item = ev_item,offset_ms = offset_ms, duration_ms = duration_ms)

                         # update 
                         subjDictList_lfo.append(E.lfo_dict)


              # convert to dataframe
              self.subjDf_lfo = pd.DataFrame(subjDictList_lfo,index=self.uElbls_thisSubjClus)
              

              # collapse data within subject (as a dict)
              self.subjDataCollapsed_lfo= {}
              self.subjDataCollapsed_lfo['psd_mean'] = np.nanmean(np.stack(self.subjDf_lfo['psd'].to_numpy()),axis=0)
              self.subjDataCollapsed_lfo['freqs'] = self.subjDf_lfo.iloc[0]['freqs']

              # store cached data in self for plotting functions
              self.lfo_psds = np.vstack(self.subjDf_lfo['psd'].to_numpy())
              self.lfo_peaks_lfoOnly = np.vstack(self.subjDf_lfo['LFO'].to_numpy())[:,0]
              self.lfo_peaks_all = np.vstack(self.subjDf_lfo['peaks'].to_numpy())[:,0]

          def plot_lfo_psd(self, ax = None):
              if ax == None:
                         f = plt.figure()
                         ax = plt.subplot(111)

              psds = self.lfo_psds
              psds_log = np.log(psds)
              psd_log_mean = np.nanmean(psds_log,axis =0)
              psd_log_sem = stats.sem(psds_log,axis =0)
              freqs = self.subjDf_lfo.iloc[0]['freqs']


              plt.plot(freqs,psd_log_mean,linewidth=3)
              plt.fill_between(freqs,psd_log_mean+psd_log_sem,psd_log_mean-psd_log_sem,alpha=0.4)
              plt.xlabel('Frequency (Hz)')
              plt.ylabel('log power')
              ax.axis('tight')



          def plot_lfo_peaks(self, lfo_only = True, ax = None,bins = 50,title=None):
              # plots distribution of lfo peaks cached fooof data. if lfo_only is True it only uses LFO peaks, if False, it plots all peaks identified (can be mutliple peaks from a single electrode) Class Subject. Overwrites Electrode.plot_lfo
              if ax == None:
                         f = plt.figure()
                         ax = plt.subplot(111)

                         if lfo_only==True:
                                  peaks = self.lfo_peaks_lfoOnly
                         else:
                                  #This uses all peaks
                                  peaks = self.lfo_peaks_all

                         plt.hist(peaks,bins = bins)
                         plt.xlabel('Frequency (Hz)')
                         plt.ylabel('count')
                         plt.title(title)

          ##############################
          #SYNCHRONY
          ##############################
          def calc_phaseSync_xTrial(self, fRange=(3,8),ev_type = None, ev_item = None,time_bin_size_ms = 100, method = 'filt_hilb',wave_numfreqs = None, wave_numcycles = None):
                # Subject class overwrites Electrode.calc_phaseSync_xTrial This function computes pairwise synchrony for all pairs in self (accounting for subject and region specifications). Creates a dataframe that can be aggregated across subjects

                # create save file name
                fname = 'xTrial-'+str(self.subj)+'-clusId-'+str(self.clusId_list)+'-ROI-'+str(self.roi_list)+'-'+str(self.paramsDict)+'-'+str(fRange)+'-'+str(ev_type)+'-'+str(ev_item)+'-'+str(time_bin_size_ms)+'-'+str(method)+'-'+str(wave_numfreqs)+'-'+str(wave_numcycles)

                # check if it is saved, then load
                fpath = scratch_dir+'/'+fname
                if (os.path.exists(fpath)==True):
                     self.subjDf_phaseSync_xTrial  = self.load_pickle(fpath)
                else:
                     #compute it and save
                     # get pair list
                     if hasattr(self,'subjPair_dict') == False:
                          self.make_pairs()

                     # init dict list
                     sync_dict_list = []
                     # loop through pairs and calculate phase sync
                     for p in np.arange(0,len(self.subjPair_dict['pair_uElbl'])):

                          # init electrode
                          E = Electrode(self.subjPair_dict['pair_uElbl_self'][p])

                          #
                          O = Electrode(self.subjPair_dict['pair_uElbl_other'][p])

                          # calc syncrony
                          E.calc_phaseSync_xTrial(other = O, fRange=fRange,ev_type = ev_type, ev_item = ev_item,time_bin_size_ms = time_bin_size_ms, method = method,wave_numfreqs = wave_numfreqs, wave_numcycles = wave_numcycles)

                          # append
                          sync_dict_list.append(E.phsSyncDict_xTrial)

                          print(p,'/',len(self.subjPair_dict['pair_uElbl']))

                     # convert to dataframe
                     self.subjDf_phaseSync_xTrial = pd.DataFrame(sync_dict_list,index = self.subjPair_dict['pair_uElbl'])     

                     #save pickle
                     self.save_pickle(self.subjDf_phaseSync_xTrial,fpath)

          def calc_phaseSync_xTime(self, fRange=(3,8),ev_type = None, ev_item = None, method = 'filt_hilb',wave_numfreqs = None, wave_numcycles = None,offset_ms_list=[0], duration_ms=1000,n_iters = 10000):
                # Subject class overwrites Electrode.calc_phaseSync_xTime This function computes pairwise synchrony for all pairs in self (accounting for subject and region specifications). Creates a dataframe that can be aggregated across subjects

                # create save file name
                fname = 'xTime-'+str(self.subj)+'-clusId-'+str(self.clusId_list)+'-ROI-'+str(self.roi_list)+'-'+str(self.paramsDict)+'-'+str(fRange)+'-'+str(ev_type)+'-'+str(ev_item)+'-'+str(method)+'-'+str(wave_numfreqs)+'-'+str(wave_numcycles)+'-'+str(offset_ms_list)+'-'+str(duration_ms)+'-'+str(n_iters)

                # keys to del (to save space)
                keys_to_del=['phsmat_clip_list', 'xvalms_clip_list', 'RVL','RVL_null','omn_p','omn_p_vsNull', 'omn_z_vsNull','RVL_zByNull']

                # check if it is saved, then load
                fpath = scratch_dir+'/'+fname
                if (os.path.exists(fpath)==True):
                     self.subjDf_phaseSync_xTime  = self.load_pickle(fpath)
                else:
                     # compute it and save
                     # get pair list
                     if hasattr(self,'subjPair_dict') == False:
                          self.make_pairs()

                     # init dict list
                     sync_dict_list = []
                     # loop through pairs and calculate phase sync
                     for p in np.arange(0,len(self.subjPair_dict['pair_uElbl'])):

                          # init electrode
                          E = Electrode(self.subjPair_dict['pair_uElbl_self'][p])

                          #
                          O = Electrode(self.subjPair_dict['pair_uElbl_other'][p])

                          # calc syncrony
                          E.calc_phaseSync_xTime(other = O, fRange=fRange,ev_type = ev_type, ev_item = ev_item, method = method,wave_numfreqs = wave_numfreqs, wave_numcycles = wave_numcycles,offset_ms_list=offset_ms_list, duration_ms=duration_ms,n_iters = n_iters)

                          # delete keys
                          for k in keys_to_del:
                                del E.phsSyncDict_xTime[k];

                          # append
                          sync_dict_list.append(E.phsSyncDict_xTime)

                          print(p,'/',len(self.subjPair_dict['pair_uElbl']))

                     # convert to dataframe
                     self.subjDf_phaseSync_xTime = pd.DataFrame(sync_dict_list,index = self.subjPair_dict['pair_uElbl'])

                     #save pickle
                     self.save_pickle(self.subjDf_phaseSync_xTime,fpath)



##############################################################################
#GROUP CLASS
###############################################################################


class Group(Subject):
# performs group-level analyses across subjects. Inherits methods and attributes from Subject class. Takes a list of subjects as an input. Has options to filter by cluster Id (when making subject object). Will have the ability to perform across electrode tests (using subjDf field) or collapsed (across subject test using subjDataCollapsed field)


          # constructor
          def __init__(self,subj_list=None,clusId_list = None,roi_list = None, paramsDict = None):
              # clusId_list... list of clusId values to include (e.g. [1,2,3] includes electrodes that were assigned to the first three clusters)
              # paramsDict .. this is optional and will overwrite default params in electrode class

              # if subj_list = None
              if subj_list==None:
              	subj_list = np.unique(groupDf['subj'])

              self.subj_list = subj_list
              self.clusId_list = clusId_list
              self.roi_list = roi_list
              self.paramsDict = paramsDict

              #initialize Subject object with first electrode so we can access methods
              Subject.__init__(self,subj = subj_list[0],clusId_list = clusId_list,roi_list=roi_list,paramsDict=paramsDict)
          def collapseBySubj_1d(self,subj_list,x):
              #This function collapses a 1d array within subject
              # subj_list ... list of subjects matching indices of array(x)

              uSubj_list = np.unique(subj_list)
              x_grp = np.zeros((len(uSubj_list)))

              for i in np.arange(0,len(uSubj_list)):
                s_idx = np.array(subj_list)==uSubj_list[i]
                x_grp[i] = np.nanmean(x[s_idx])

              return x_grp

          


          ##############################
          #ERP
          ##############################

          def calc_erp(self,ev_type = None, ev_item = None,time_bin_size_ms = 100,apply_gauss_smoothing = True,gauss_sd_scaling = 0.1):
              # loops through electrodes and updates electrode dict. Group class; overwrites S.calc_erp

              # initialize empty dataframe to collect electrode..-level data
              self.groupDf_erp = pd.DataFrame([])
              groupDataCollapsed_list = []

              #loop through electrode list
              for s in self.subj_list:
                         S = Subject(subj=s,clusId_list = self.clusId_list,roi_list=self.roi_list,paramsDict = self.paramsDict)

                         #if no matching electrodes, continue
                         if S.isempty == True:
                                  continue

                         # calcERP
                         S.calc_erp(ev_type = ev_type, ev_item = ev_item,time_bin_size_ms = time_bin_size_ms,apply_gauss_smoothing = apply_gauss_smoothing,gauss_sd_scaling = gauss_sd_scaling)

                         # update groupDf_erp
                         self.groupDf_erp=self.groupDf_erp.append(S.subjDf_erp)

                         # update 
                         groupDataCollapsed_list.append(S.subjDataCollapsed_erp)


              # convert collapsed data to data frame
              self.groupDfCollapsed_erp = pd.DataFrame(groupDataCollapsed_list,index=np.unique(self.groupDf_erp['subj']))

          def plot_erp(self,ax = None,figsize=None,fsize_tick=12,alpha=0.5,xRange_ms = (-500,1000),use_collapsed = False,yL=None):
              # for class group. overwrites Electrode.plot_erp and Subject.plot_erp
              # It only plots z-scored data that has been time binned
              if ax == None:
                         f = plt.figure(figsize=figsize)              
                         ax = plt.subplot(111)
              # get xval, mean and sem
              xval = self.groupDf_erp['erp_xval_ms'][0]
              if use_collapsed == False:
                         # use electrode-level data
                         erp_mat = np.stack(self.groupDf_erp['abs_z_erp_mean_ms'].to_numpy()) 
                         erp_mean = np.nanmean(erp_mat,axis=0)
                         erp_sem = stats.sem(erp_mat,axis=0)
              else:
                         # use collapsed data
                         erp_mat = np.stack(self.groupDfCollapsed_erp['abs_z_erp_mean_ms'].to_numpy()) 
                         erp_mean = np.nanmean(erp_mat,axis=0)
                         erp_sem = stats.sem(erp_mat,axis=0)


              # plot it
              plt.plot(xval,erp_mean,alpha=alpha,label = str(self.groupDfCollapsed_erp['erp_params'][0]['ev_type'])+'/'+str(self.groupDfCollapsed_erp['erp_params'][0]['ev_item']))
              plt.fill_between(xval,erp_mean+erp_sem,erp_mean-erp_sem,alpha=0.5)

              # use time bins data
              xt = np.array([xRange_ms[0],0,xval[np.argmin(np.abs(xval-500))],xRange_ms[1]])
              ax.set_xticks(xt)
              ax.set_xticklabels(xt,fontsize=fsize_tick)
              ax.set_xlim(xRange_ms[0],xRange_ms[1])

              if yL!=None:
                         ax.set_ylim(yL[0],yL[1])



          def plot_erp_byType(self,ev_item=None,use_collapsed=False,apply_gauss_smoothing = True,gauss_sd_scaling = 0.1,yL=None):
              # class Group: overwrites Electrode.plot_erp_byType
              # if ev_item = None, it will combine data from all items. Otherwise, it will only plot data from the specified item (left/right/mouth)
              f = plt.figure()
              ax = plt.gca()
              evType_list = np.unique(groupDf['ev_type'][0])

              for t in evType_list:
                         self.calc_erp(ev_type=t,ev_item = ev_item,apply_gauss_smoothing = apply_gauss_smoothing,gauss_sd_scaling = gauss_sd_scaling)
                         self.plot_erp(ax=ax,use_collapsed=use_collapsed,yL=yL)

              plt.legend()
              plt.vlines(x=0,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],alpha = 0.3,color='r')
              plt.xlabel('Time from stimulus (ms)')
              plt.ylabel('absolute value of z-scored erp')
              plt.title('Clusters = '+str(self.clusId_list)+', Regions = '+str(self.roi_list)+'; N elec(subj) = '+str(len(self.groupDf_erp))+' ('+str(len(self.groupDfCollapsed_erp))+')')          

          def plot_erp_byItem(self,ev_type = 'MOVE',use_collapsed=False,apply_gauss_smoothing = True,gauss_sd_scaling = 0.1,yL=None):
              # class Group: overwrites Electrode.plot_erp_byItem

              f = plt.figure()
              ax = plt.gca()

              evItem_list = np.unique(groupDf['ev_item'][0])

              for i in evItem_list:
                         self.calc_erp(ev_type=ev_type,ev_item=i,apply_gauss_smoothing = apply_gauss_smoothing,gauss_sd_scaling = gauss_sd_scaling)
                         self.plot_erp(ax=ax,use_collapsed=use_collapsed,yL=yL)                                  
              plt.legend()
              plt.vlines(x=0,ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1],alpha = 0.3,color='r')
              plt.xlabel('Time from stimulus (ms)')
              plt.ylabel('absolute value of z-scored erp')
              plt.title('Clusters = '+str(self.clusId_list)+', Regions = '+str(self.roi_list)+'; N elec(subj) = '+str(len(self.groupDf_erp))+' ('+str(len(self.groupDfCollapsed_erp))+')')          


          ##############################
          #LFO
          ##############################

          def calc_lfo(self,ev_type = None, ev_item = None,offset_ms = 0, duration_ms = 1000):
              # loops through subjects and updates groupDf. Group class, overwrites Subject.calc_lfo


              # initialize empty dataframe to collect electrode..-level data
              self.groupDf_lfo = pd.DataFrame([])
              groupDataCollapsed_list = []


              #loop through electrode list
              for s in self.subj_list:
                         S = Subject(subj=s,clusId_list = self.clusId_list,roi_list=self.roi_list,paramsDict = self.paramsDict)

                         #if no matching electrodes, continue
                         if S.isempty == True:
                                  continue

                         # calcERP
                         S.calc_lfo(ev_type = ev_type, ev_item = ev_item,offset_ms = offset_ms, duration_ms = duration_ms)

                         # update groupDf_erp
                         self.groupDf_lfo = self.groupDf_lfo.append(S.subjDf_lfo)

                         # update 
                         groupDataCollapsed_list.append(S.subjDataCollapsed_lfo)


              # convert collapsed data to data frame
              self.groupDfCollapsed_lfo = pd.DataFrame(groupDataCollapsed_list,index=np.unique(self.groupDf_lfo['subj']))


              # store cached data in self for plotting functions
              #self.lfo_psds = np.vstack(self.groupDf_lfo['psd'].to_numpy())
              self.lfo_peaks_lfoOnly = np.vstack(self.groupDf_lfo['LFO'].to_numpy())[:,0]
              self.lfo_peaks_all = np.vstack(self.groupDf_lfo['peaks'].to_numpy())[:,0]

          def plot_lfo_peaks(self, lfo_only = True, ax = None,bins = 50,title = None):
              # overwrites S.plot_lfo_peaks. Adds a fig title
              if ax == None:
                         f = plt.figure()
                         ax = plt.subplot(111)

              if lfo_only==True:
                peaks = self.lfo_peaks_lfoOnly
              else:
                #This uses all peaks
                peaks = self.lfo_peaks_all


              plt.hist(peaks,bins = bins,color='0.5')
              plt.xlabel('Frequency (Hz)')
              plt.ylabel('count')
              if title==None:
                      plt.title('Clusters = '+str(self.clusId_list)+', Regions = '+str(self.roi_list)+'; N elec(subj) = '+str(len(self.groupDf_lfo))+' ('+str(len(self.groupDfCollapsed_lfo))+')')
              else:
                      plt.title(title)
          ##############################
          #SYNCHRONY
          ##############################
          def calc_phaseSync_xTrial(self, fRange=(3,8),ev_type = None, ev_item = None,time_bin_size_ms = 100, method = 'filt_hilb',wave_numfreqs = None, wave_numcycles = None):
                # Group class overwrites Subject.calc_phaseSync_xTrial This function computes pairwise synchrony for all pairs in self (accounting for subject and region specifications). Creates a dataframe that can be aggregated across subjects


                # init df
                self.groupDf_phaseSync_xTrial = pd.DataFrame()
                # loop through subjects and calculate phase sync
                for s in self.subj_list:

                     # Init subj
                     S = Subject(subj=s,clusId_list = self.clusId_list, roi_list = self.roi_list, paramsDict = self.paramsDict)

                     # calc syncrony
                     S.calc_phaseSync_xTrial(fRange=fRange,ev_type = ev_type, ev_item = ev_item,time_bin_size_ms = time_bin_size_ms, method = method,wave_numfreqs = wave_numfreqs, wave_numcycles = wave_numcycles)

                     # append
                     self.groupDf_phaseSync_xTrial=self.groupDf_phaseSync_xTrial.append(S.subjDf_phaseSync_xTrial)


          def calc_phaseSync_xTime(self, fRange=(3,8),ev_type = None, ev_item = None, method = 'filt_hilb',wave_numfreqs = None, wave_numcycles = None,offset_ms_list=[0], duration_ms=1000,n_iters = 10000):
                #  Group class overwrites Subject.calc_phaseSync_xTime This function computes pairwise synchrony for all pairs in self (accounting for subject and region specifications). Creates a dataframe that can be aggregated across subjects

                # init Df
                self.groupDf_phaseSync_xTime = pd.DataFrame()

                # loop through subjects and calculate phase sync
                for s in self.subj_list:

                     # init electrode
                     S = Subject(subj=s,clusId_list = self.clusId_list, roi_list = self.roi_list, paramsDict = self.paramsDict)

                     # calc syncrony
                     S.calc_phaseSync_xTime(fRange=fRange,ev_type = ev_type, ev_item = ev_item, method = method,wave_numfreqs = wave_numfreqs, wave_numcycles = wave_numcycles,offset_ms_list=offset_ms_list, duration_ms=duration_ms,n_iters = n_iters)

                     # append
                     self.groupDf_phaseSync_xTime=self.groupDf_phaseSync_xTime.append(S.subjDf_phaseSync_xTime)

                # Add theta columns (for convenience)
                # has_theta_self
                has_theta_self = np.isnan(self.groupDf_phaseSync_xTime['theta_cf'].to_numpy())==False
                self.groupDf_phaseSync_xTime.insert(loc=len(self.groupDf_phaseSync_xTime.columns),column='has_theta_self',value=has_theta_self) 

                # has_theta_other
                has_theta_other = np.isnan(self.groupDf_phaseSync_xTime['other_theta_cf'].to_numpy())==False
                self.groupDf_phaseSync_xTime.insert(loc=len(self.groupDf_phaseSync_xTime.columns),column='has_theta_other',value=has_theta_other) 

                # theta_cf_diff
                theta_cf_diff = np.absolute((np.squeeze(self.groupDf_phaseSync_xTime['theta_cf'].to_numpy()))-(np.squeeze(self.groupDf_phaseSync_xTime['other_theta_cf'].to_numpy())))
                self.groupDf_phaseSync_xTime.insert(loc=len(self.groupDf_phaseSync_xTime.columns),column='theta_cf_diff',value=theta_cf_diff) 
          def queryGroupPhsDf(self,phsdf, query_str = None):
            # This function returns a  copy of groupDf_phaseSync_xTime that has been queried by the query str. Can use the returned copy in additional functions below
            # inputs
            #phsdf ... group phase data frame (e.g., a copy of groupDf_phaseSync_xTime). Can pass in queried da

            if query_str!=None:
              phsdf =phsdf.query(query_str)

            return phsdf

          def updateClusDataInGroupPhsDf(self,phsdf):
            # This function updates the clustering data in a group phase data frame (e.g., a copy of groupDf_phaseSync_xTime) with the current matlab-generated clustering data (clusId, postMove_HFA, etc). This allows you to update the matlab clustering algorithm without re-running pairwise phase data
            # inputs
            #phsdf ... group phase data frame (e.g., a copy of groupDf_phaseSync_xTime). Can pass in queried data

            # returns updated phsdf data

            # identify columns to update (automatically updates self and other_ values for each pair)
            col_to_update = ['clusId','postInst_HFA', 'postMove_HFA','postInst_theta','effectorSelective', 'instructionSelective','postMove_theta','postInst_beta','postMove_beta','postInst_LFA','postMove_LFA']

            # loop through each pair
            for idx,row in phsdf.iterrows():

                # get E data
                E_uElbl = row['uElbl']
                E_dict = groupDf.query('uElbl == @E_uElbl').iloc[0]   

                # get other_uElbl data
                O_uElbl = row['other_uElbl']
                O_dict = groupDf.query('uElbl == @O_uElbl').iloc[0]

                # loop through col to update
                for c in col_to_update:

                    # skip if the column is not included in phsdf
                    if (c in phsdf.columns) == False:
                        continue
                    # update phsdf with electrode data
                    phsdf.at[idx,c] = E_dict[c]

                    # update phsdf with other_ data
                    phsdf.at[idx,'other_'+c] = O_dict[c]

            return phsdf
          def getVal_phaseSync_xTime(self,lbl, offset_of_interest_ms,phsdf=None):
            # This function returns the value of a particular field across all pairs in phsdf for the offset of interest
            if phsdf is None:
               phsdf = self.groupDf_phaseSync_xTime.copy()


            # get time interval of interest
            t_idx = np.argmin(np.abs(phsdf.iloc[0]['params']['offset_ms'] - offset_of_interest_ms))

            # get data
            a = np.vstack(phsdf[lbl].to_numpy())
            val = a[:,t_idx]


            return val
           



          def plot_phaseSync_xTime_rvl_tstat_ByTime(self,phsdf = None,ax = None,offset_of_interest_ms=None):
            # Plot RVL over time
            # Inputs
            #phsdf ... option to provide a groupDf that has been queried. If None will revert to cached data
            #offset_of_interest_ms ... option to mark the time window we are interested in
            
            if ax == None:
                f = plt.figure()
                ax = plt.subplot(111)
            if phsdf is None:
                phsdf = self.groupDf_phaseSync_xTime.copy()

            x = phsdf.iloc[0]['xval_ms']
            a = np.vstack(phsdf['RVL_tVsNull_t'].to_numpy())
            a_m = np.nanmean(a,axis=0)
            a_sem = stats.sem(a,axis=0,nan_policy='omit')
            ax.plot(x,a_m,linewidth=3,alpha=0.5,color='C0')
            ax.fill_between(x,a_m+a_sem, a_m-a_sem,alpha=0.5)
            ax.axis('tight')
            ax.set_ylim((0,ax.get_ylim()[1]))
            ax.vlines(0,ax.get_ylim()[0],ax.get_ylim()[1],linestyle='dashed',color='0.5',alpha=0.5)
            ax.set_xlabel('Time from cue (ms)')
            ax.set_ylabel('t-statistic')
            ax.set_title('n = '+str(len(phsdf)))
            if offset_of_interest_ms!=None:
                t_idx = np.argmin(np.abs(phsdf.iloc[0]['params']['offset_ms'] - offset_of_interest_ms))
                ax.vlines(x[t_idx],ax.get_ylim()[0],ax.get_ylim()[1],linestyle='dashed',color = 'C1')
          def plot_phaseSync_xTime_tstat_dist(self,phsdf = None,ax= None,offset_of_interest_ms = None):
            # plot RVL_vs_null t-stat distribution
            # Inputs
            #phsdf ... option to provide a groupDf that has been queried. If None will revert to cached data
            #offset_of_interest_ms ... specify the time window we are interested in
            if ax == None:
                f = plt.figure()
                ax = plt.subplot(111)
            if phsdf is None:
                phsdf = self.groupDf_phaseSync_xTime.copy()

            t_idx = np.argmin(np.absolute(phsdf['params'].iloc[0]['offset_ms']-offset_of_interest_ms))
            rvl_tstats = np.stack(phsdf['RVL_tVsNull_t'].to_numpy())
            rvl_tstats_thisOffset = rvl_tstats[:,t_idx]
            ax.hist(rvl_tstats_thisOffset,color='0.5',bins=100)
            t,p =stats.ttest_1samp(rvl_tstats_thisOffset,0)
            ax.set_title('t = '+str(np.round(t,2))+' p = '+str(np.round(p,4)))
            ax.set_xlabel('t statistic')
            ax.set_ylabel('number of electrode pairs')

            return rvl_tstats_thisOffset




