#Ke Chen
import obspy
import pyflex

obs_st = obspy.read("../pyflex/data/BK.WENL.BHE.sac")   # obs_st is a stream object
syn_st = obspy.read("../pyflex/data/BK.BH.WENL.e")      # syn_st is a stream object
syn_st.resample(obs_st[0].stats.sampling_rate)



#obs_st.trim(starttime=max(obs_st[0].stats.starttime, syn_st[0].stats.starttime), endtime=min(obs_st[0].stats.endtime, syn_st[0].stats.endtime))
#syn_st.trim(starttime=max(obs_st[0].stats.starttime, syn_st[0].stats.starttime), endtime=min(obs_st[0].stats.endtime, syn_st[0].stats.endtime))


# preprocessing
obs_st.detrend("linear")
obs_st.taper(max_percentage=0.05, type="hann")
obs_st.filter("lowpass",freq=0.25,corners=2,zerophase=True)

syn_st.detrend("linear")
syn_st.taper(max_percentage=0.05, type="hann")
syn_st.filter("lowpass",freq=0.25,corners=2,zerophase=True)

config = pyflex.Config(
                       min_period=4.0, max_period=100.0,
                       stalta_waterlevel=0.08, tshift_acceptance_level=15.0,
                       dlna_acceptance_level=1.0, cc_acceptance_level=0.80,
                       c_0=0.7, c_1=4.0, c_2=0.0, c_3a=1.0, c_3b=2.0, c_4a=3.0, c_4b=10.0)

windows = pyflex.select_windows(obs_st, syn_st, config, plot=True)


