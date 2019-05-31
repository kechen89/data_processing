# 1_event_catalog_waveform_search.py
"""
Seismic event catalog request for NCEDC

cat_ncss has the standard processing ('reviewed' and 'final') origins from the NCSS catalog
    I also collect the focal mechanisms in this catalog
    
cat_dd has the Double-Difference origins from the DD catalog

"""

import os, sys
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.imaging.beachball import beachball
import numpy as np

######################################
# Main Program

# Event search parameters 
# Time period limits
tstart = UTCDateTime("2000-01-01")
tend = UTCDateTime("2018-12-31")

# Geographic limits
min_lon = -123.10
max_lon = -121.20
min_lat = 37.20
max_lat = 38.60

# Magnitude limits
min_mag = 4.0
max_mag = 5.0

# Data center
data_center = 'NCEDC'

# SW4 parameters
#   t0 = time shift before origin time, 
#       simulation starts this many seconds before the origin time
#   duration = source duration, should be short for impulsive source-time function
#   freq = frequancy parameter of source (radians/second), based on duration
#   seismogram_time = maximum time of the simulation
t0 = 5.0
duration = 1.0
freq = 2*3.14159265359/duration
seismogram_time = 60.

######################################

client = Client(data_center)                 # initialize a client object
print ('... running query ...')
cat_ncss = client.get_events(starttime=tstart, endtime=tend,             # retrieve event data from sever, return a Catlog object
                              minlatitude=min_lat, minlongitude=min_lon, 
                              maxlatitude=max_lat, maxlongitude=max_lon,
                              minmagnitude=min_mag, maxmagnitude=max_mag,
                              catalog='NCSS', includemechanisms=True)

cat_dd = client.get_events(starttime=tstart, endtime=tend,
                              minlatitude=min_lat, minlongitude=min_lon, 
                              maxlatitude=max_lat, maxlongitude=max_lon,
                              minmagnitude=min_mag, maxmagnitude=max_mag,
                              catalog='DD')

# Report number of events
print ('cat_ncss found: '+str(cat_ncss.count())+' event(s)')   
print ('cat_dd found:   '+str(cat_dd.count())+' event(s)')  

# Plot NCSS events
fig = cat_ncss.plot(method="basemap", label='magnitude', projection='local', color='date', resolution='i', title='NCSS Event Catalog') 
fig.savefig('catalog_ncss_map.png', figsize=(6,6), dpi=500)

# Write out catalogs for all events
filename = 'catalog_ncss.xml'
cat_ncss.write(filename, format="QUAKEML")
filename = 'catalog_dd.xml'
cat_dd.write(filename, format="QUAKEML")

catalog_file = './catalog_ncss.txt'
if os.path.exists(catalog_file)==True:
    os.remove(catalog_file)
f=open(catalog_file,'w')
f.write('# lon, lat, depth_km, mag, time, event_ID\n')
for event in cat_ncss.events:
    #event=cat_ncss.events[k]
    evla=event.origins[0].latitude
    evlo=event.origins[0].longitude
    depth=event.origins[0].depth/1000
    mag=event.magnitudes[0].mag
    time=str(event.origins[0].time)
    event_id=str(event.resource_id).split('/')[-1]
    line='%.4f\t%.4f\t%8.2f\t%.2f\t%s\t%s\n' % (evlo,evla,depth,mag,time,event_id)
    f.write(line)
f.close()

f_mt_psmeca = open('cat_mt_psmeca.txt', 'w')
f_dc_psmeca = open('cat_dc_psmeca.txt', 'w')

# Now loop over events in the cat_ncss
for event in cat_ncss.events:
    event_id = str(event.resource_id).split('/')[-1]
    print('#################')
    print('Event_id: ', event_id)
    if not os.path.isdir(event_id):
        os.mkdir(event_id)
    f = open(event_id+'/'+event_id+'.event_report', 'w')

    # Find NCSS preferred origin (evaluation_mode = 'final'), save as dictionary
    ncss_origin = {}
    for origin in event.origins:
        if origin.evaluation_status == 'final':
        #if origin.evaluation_status == 'reviewed':
            ncss_origin['evid'] = event_id
            ncss_origin['latitude'] = origin.latitude
            ncss_origin['longitude'] = origin.longitude
            ncss_origin['depth'] = origin.depth
            ncss_origin['otime'] = str(origin.time)
            #ncss_df = pd.DataFrame.from_dict(ncss_origin)
            f_ncss = open(event_id+'/'+str(event_id)+'.ncssorigin', 'w')
            f_ncss.write(str(ncss_origin))
            f_ncss.close()
            f.writelines('NCSS origin: \n')
            f.writelines(str(ncss_origin))
            f.writelines('\n')
    # Find DD origin, save as dictionary, write to file
    # Note: this assumes there is only one DD origin for this event
    dd_origin = {}
    for event_dd in cat_dd.events:
        event_id_dd = str(event_dd.resource_id).split('/')[-1]
        if event_id_dd == event_id:
            print('Event_id_dd: ', event_id_dd)
            #print('Event_dd_origin: ', event_dd.origins, len(cat_dd.events[0].origins))
            dd_origin['evid'] = event_id
            dd_origin['latitude'] = event_dd.origins[0].latitude 
            dd_origin['longitude'] = event_dd.origins[0].longitude 
            dd_origin['depth'] = event_dd.origins[0].depth 
            dd_origin['otime'] = str(event_dd.origins[0].time)
            f_dd = open(event_id+'/'+str(event_id)+'.ddorigin', 'w')
            f_dd.write(str(dd_origin))
            f_dd.close()
            f.writelines('DD origin: \n')
            f.writelines(str(dd_origin))
            f.writelines('\n')
            t_sw4time = open(event_id+'/'+str(event_id)+'.sw4time', 'w')
            t_sw4time.writelines('# DD origin time: '+dd_origin['otime']+'\n')
            t_sw4time.writelines('time t='+str(seismogram_time)+' utctime='+(event_dd.origins[0].time - t0).strftime('%d-%m-%Y:%H:%M:%S.%f'))
            t_sw4time.close()
            break
    # Find preferred magnitude
    for magnitude in event.magnitudes:
        if magnitude.resource_id == event.preferred_magnitude_id:
            mag = magnitude.mag
            m0 = ( 10.0**((3./2.)*(mag + 10.7)) ) / 1e7
            f.writelines('preferred mag: '+str(round(mag,2))+'\n')
            
    # Find distance between NCSS final and DD locations
    if ncss_origin != {} and dd_origin != {}:
        epi_dist_m = gps2dist_azimuth(lat1=ncss_origin['latitude'], lon1=ncss_origin['longitude'], lat2=dd_origin['latitude'], lon2=dd_origin['longitude'])
        print('Distance between NCSS final and DD epicenters: ', epi_dist_m[0] )
        f.write('Distance (meters) between NCSS final and DD epicenters: '+str(round(epi_dist_m[0],2)))
        f.writelines('\n')
    
    # Find focal mechanisms, check for double-couple (dc) and moment tensor (mt)
    for focal_mechanism in event.focal_mechanisms:
        ft_type = None
        if focal_mechanism.nodal_planes != None:
            fm_type = 'dc'
            fm_dc = {
                    'strike': focal_mechanism.nodal_planes.nodal_plane_1.strike,
                    'dip': focal_mechanism.nodal_planes.nodal_plane_1.dip,
                    'rake': focal_mechanism.nodal_planes.nodal_plane_1.rake,
                    'mag': mag,
                    'm0': m0
                    }
            dc_arr = np.array( [ fm_dc['strike'], fm_dc['dip'], fm_dc['rake'] ])
            dc = dc_arr.tolist()
            beachball(dc, size=200, linewidth=2, facecolor='b')
            line = 'source lat='+str(round(dd_origin['latitude'],5))
            line += ' lon='+str(round(dd_origin['longitude'],5))
            line += ' depth='+str(round(dd_origin['depth'],1))
            line += ' t0='+str(round(t0,1))
            line += ' m0='+str(round(fm_dc['m0'],1))
            line += ' strike='+str(round(fm_dc['strike'],0))
            line += ' dip='+str(round(fm_dc['dip'],0))
            line += ' rake='+str(round(fm_dc['rake'],0))
            line += ' type=Gaussian freq='+str(round(freq,3))
            f_sw4source = open(event_id+'/'+event_id+'.sw4source_dc', 'w')
            f_sw4source.writelines(line)
            f_sw4source.close()
            f_psmeca_dc = open(event_id+'/'+event_id+'.psmeca_dc', 'w')
            line = str(round(dd_origin['longitude'],5))+' '+str(round(dd_origin['latitude'],5))+' '
            line += str(round(dd_origin['depth']/1000.,1))+' '+str(round(fm_dc['strike'],0))+' '
            line += str(round(fm_dc['dip'],0))+' '+str(round(fm_dc['rake'],0))+' '
            line += str(round(fm_dc['mag'],2))+' '
            line += str(round(dd_origin['longitude']-0.1,5))+' '+str(round(dd_origin['latitude'],5))+' '
            line += event_id+'\n'
            f_psmeca_dc.writelines(line)
            f_psmeca_dc.close()
            f.writelines('Double-couple focal mechanism: \n')
            f.writelines(str(fm_dc))
            f.writelines('\n')
            
        if focal_mechanism.moment_tensor != None:
            fm_type = 'mt'
            fm_mt = {
                    'm_rr': focal_mechanism.moment_tensor.tensor['m_rr']/focal_mechanism.moment_tensor.scalar_moment,
                    'm_tt': focal_mechanism.moment_tensor.tensor['m_tt']/focal_mechanism.moment_tensor.scalar_moment,
                    'm_pp': focal_mechanism.moment_tensor.tensor['m_pp']/focal_mechanism.moment_tensor.scalar_moment,
                    'm_rt': focal_mechanism.moment_tensor.tensor['m_rt']/focal_mechanism.moment_tensor.scalar_moment,
                    'm_rp': focal_mechanism.moment_tensor.tensor['m_rp']/focal_mechanism.moment_tensor.scalar_moment,
                    'm_tp': focal_mechanism.moment_tensor.tensor['m_tp']/focal_mechanism.moment_tensor.scalar_moment,
                    'scalar_moment': focal_mechanism.moment_tensor.scalar_moment
                    }
            fm_mt['m0'] = fm_mt['scalar_moment']
            fm_mt['mag'] = (2/3) * np.log10(fm_mt['m0']*1.e7) - 10.73

            mt_arr = np.array([fm_mt['m_rr'], fm_mt['m_tt'], fm_mt['m_pp'], fm_mt['m_rt'], fm_mt['m_rp'], fm_mt['m_tp']])*fm_mt['m0']*1.e7
            mt = mt_arr.tolist()
            beachball(mt, size=200, linewidth=2, facecolor='r')
            line = 'source lat='+str(round(dd_origin['latitude'],5))
            line += ' lon='+str(round(dd_origin['longitude'],5))
            line += ' depth='+str(round(dd_origin['depth'],1))
            line += ' t0='+str(round(t0,1))
            line += ' m0='+str(round(fm_mt['m0'],1))
            line += ' mxx='+str(round(fm_mt['m_tt'],5))+' myy='+str(round(fm_mt['m_pp'],5))
            line += ' mzz='+str(round(fm_mt['m_rr'],5))+' mxy='+str(round(-fm_mt['m_tp'],5))
            line += ' mxz='+str(round(fm_mt['m_rp'],5))+' myz='+str(round(-fm_mt['m_rp'],5))
            line += ' type=Gaussian freq='+str(round(freq,3))
            f_sw4source = open(event_id+'/'+event_id+'.sw4source_mt', 'w')
            f_sw4source.writelines(line)
            f_sw4source.close()
            f_psmeca_mt = open(event_id+'/'+event_id+'.psmeca_mt', 'w')
            line = str(round(dd_origin['longitude'],5))+' '+str(round(dd_origin['latitude'],5))+' '
            line += str(round(dd_origin['depth']/1000.,1))+' '
            line += str(round(focal_mechanism.moment_tensor.tensor['m_rr']/1e15,6))+' '
            line += str(round(focal_mechanism.moment_tensor.tensor['m_tt']/1e15,6))+' '
            line += str(round(focal_mechanism.moment_tensor.tensor['m_pp']/1e15,6))+' '
            line += str(round(focal_mechanism.moment_tensor.tensor['m_rt']/1e15,6))+' '
            line += str(round(focal_mechanism.moment_tensor.tensor['m_rp']/1e15,6))+' '
            line += str(round(focal_mechanism.moment_tensor.tensor['m_tp']/1e15,6))+' '
            line += '22 '
            line += str(round(dd_origin['longitude']+0.1,5))+' '+str(round(dd_origin['latitude'],5))+' '
            line += event_id+'\n'
            f_psmeca_mt.writelines(line)
            f_psmeca_mt.close()
            f.writelines('Moment tensor focal mechanism: \n')
            f.writelines(str(fm_mt))
            f.writelines('\n')
        print ('event_id, fm_type: ', event_id, fm_type)
    f_origins = open(event_id+'/'+event_id+'.origins', 'w')
    line = 'lon lat depth mag author\n'
    if ncss_origin != {}:
        line += str(round(ncss_origin['longitude'],5))+' '+str(round(ncss_origin['latitude'],5))+' '
        line += str(round(ncss_origin['depth']/1000.,5))+' '+str(event.magnitudes[0].mag)+' NCSS'+'\n'
    if dd_origin != {}:
        line += str(round(dd_origin['longitude'],5))+' '+str(round(dd_origin['latitude'],5))+' '
        line += str(round(dd_origin['depth']/1000.,5))+' '+str(event.magnitudes[0].mag)+' DD'+'\n'
    f_origins.writelines(line)
    f_origins.close()
    # Close event_report file
    f.close()
