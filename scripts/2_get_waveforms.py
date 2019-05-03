# 2_get_waveforms.py

from obspy import read_events
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.mass_downloader import RectangularDomain, Restrictions, MassDownloader
import os
import ast

# Parameters for weaveform search
#maxradius = 1.0# not used
min_length = 0.5
tbefore = 20.0
tafter = 60.0

# Geographic limits 
min_lon = -123.10
max_lon = -121.20
min_lat = 37.20
max_lat = 38.60

#min_lon = -122.60
#max_lon = -121.80
#min_lat = 37.60
#max_lat = 38.00

# Data center
data_center = 'NCEDC'
client = Client(data_center)

# Networks and channels
#networks = ['BK', 'NC', 'NP', 'CE', 'YU']
#channels = ['BH?', 'HN?', 'HH?', 'BL?', 'CN?']
#locations = ['', '--', '00', '01', '04']

#net_chans = ['BK.BH', 'BK.HN']
#net_chans = ['BK.HH', 'BK.BL', 'BK.BN', 'BK.CL', 'BK.CN']
#net_chans = ['WR.HN', 'PG.HN', 'GS.HN', 'GM.HN']
#net_chans = ['NC.HN', 'NP.HN', 'CE.HN']
# 'GS.HN', 'GM.HN', 'WR.HN', 'PG.HN'

net_chans = ['BK.BH'] #, 'BK.HN', 'BK.BN','BK.BL', 'BK.CL', 'NC.HN', 'NP.HN', 'CE.HN']


cat_ncss = read_events('catalog_ncss.xml')

for event in cat_ncss:
    event_id = str(event.resource_id).split('/')[-1]
    print('#################')
    print('Event_id: ', event_id)
    raw_wf_data_dir = event_id+'/RAW/'
    if not os.path.isdir(raw_wf_data_dir):
        os.mkdir(raw_wf_data_dir)
    sta_metadata_dir = raw_wf_data_dir+'/_station_xml/'
    if not os.path.isdir(sta_metadata_dir):
        os.mkdir(sta_metadata_dir)
    
    # Read NCSS preferred origin (evaluation_mode = 'final') file, was saved as dictionary
    f_ncssorigin = open(event_id+'/'+event_id+'.ncssorigin', 'r')
    ncss_origin = ast.literal_eval(f_ncssorigin.readlines()[0])
        
    print ('tbefore: ', tbefore, ' tafter: ', tafter, ' tlength: ', tafter+tbefore )
    starttime=UTCDateTime(ncss_origin['otime']) - tbefore
    endtime=UTCDateTime(ncss_origin['otime']) + tafter
    print('starttime: ', starttime, ' endtime: ', endtime)
   # Set-up domain
    domain = RectangularDomain(minlatitude=min_lat, maxlatitude=max_lat,
                           minlongitude=min_lon, maxlongitude=max_lon)
#    domain = CircularDomain(latitude=ncss_origin['latitude'], longitude=ncss_origin['longitude'],
#                            minradius=0.0, maxradius=maxradius)
    
    for net_chan in net_chans:
        raw_wf_data_dir = event_id+'/RAW/'+net_chan
        if not os.path.isdir(raw_wf_data_dir):
            os.mkdir(raw_wf_data_dir)
        sta_metadata_dir = raw_wf_data_dir+'/_station_xml/'
        if not os.path.isdir(sta_metadata_dir):
            os.mkdir(sta_metadata_dir)
        net = net_chan.split('.')[0]
        chan = net_chan.split('.')[1]
        print('... working on ... ', net_chan, net, chan)
        # Here I try to fine-tune channel and location priorities
        if net_chan in [ 'BK.BH', 'BK.HH' ]:
            channels = [ chan+'[ZNE]' ]
            locations = ['00', '', '01']
#        if net_chan == 'BK.BL':
#            channels = [ chan+'[123]' ]
#            locations = [ '', '00', '--',]
        if net_chan == 'BK.BN':
            channels = [ chan+'[123]' ]
            locations = ['40', '', '00']
        if net_chan == 'BK.CN':
            channels = [ chan+'[123]' ]
            locations = ['40', '', '00']
        if net_chan == 'BK.CL':
            channels = [ chan+'[123]' ]
            locations = ['--', '', '00']
        if net_chan == 'BK.HN':
            channels = [ chan+'[ZNE]', chan+'[123]' ]
            locations = ['', '01', '00', '40']
        if net_chan in [ 'NC.HN', 'NP.HN']:
            channels = [ chan+'[ZNE]', chan+'[Z12]']
            locations = ['', '00', '01']
        if net_chan in [ 'CE.HN']:
            channels = [ chan+'[ZNE]' ]
            locations = ['10', '', '--']
        if net_chan in [ 'PG.HN' ]:
            channels = [ chan+'[ZNE]' ]
            locations = ['--']
        if net_chan in [ 'WR.HN' ]:
            channels = [ chan+'[ZNE]' ]
            locations = ['']
        if net_chan in [ 'GS.HN' ]:
            channels = [ chan+'[Z12]' ]
            locations = ['', '00', '20']
        if net_chan in [ 'GM.HN' ]:
            channels = [ chan+'[ZNE]' ]
            locations = ['01']

        # Restrictions for waveform query, this time use minimum length restriction
        restrictions = Restrictions(starttime=UTCDateTime(ncss_origin['otime']) - tbefore,
                endtime=UTCDateTime(ncss_origin['otime']) + tafter,
                network=net,
                channel_priorities=channels,
                location_priorities=locations,
                reject_channels_with_gaps=True,
                minimum_length=min_length,
                minimum_interstation_distance_in_m=1)        
        # Use massdownloader
        mdl = MassDownloader(providers=[data_center])
        
        try:
            mdl.download(domain, restrictions, mseed_storage=raw_wf_data_dir,
                         stationxml_storage=sta_metadata_dir)
        except:
            print ('\n\n#########     DOWNLOAD ERROR: moving on      ##########\n\n')

    
