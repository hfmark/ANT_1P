import numpy as np
from obspy import read_inventory, UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core.inventory.util import Azimuth, SampleRate
from obspy.io.xseed.core import _read_seed
from glob import glob
from copy import deepcopy
import os, sys

####
# make a single stationxml(->dataless) file for all networks and instruments
# get IRIS stations from IRIS, and ENAP from disk
####

inv_all = Inventory(networks=[],source='IRIS_ENAP')

client = Client('IRIS')

# C, C1, G, XB
inv = client.get_stations(network='C,C1,G,XB', channel='LH?', level='response',\
            minlatitude=-55.5, maxlatitude=-43.1,\
            minlongitude=-76.3, maxlongitude=-65.3)

for net in inv.networks:
    inv_all.networks.append(net)  # add all of the simpler networks


# 1P, which requires GUMN and OHRS be manually set to T120
inv = client.get_stations(network='1P', channel='LH?', level='response',\
            minlatitude=-55.5, maxlatitude=-43.1,\
            minlongitude=-76.3, maxlongitude=-65.3)

horizon_sta = ['GUMN','OHRS']
GUA = inv.networks[0].select(station='X').copy()
for ista in inv.networks[0].stations:
    sta = inv.networks[0].select(station=ista.code).stations[0]
    if ista.code in horizon_sta:
        newsta = _read_seed('seed/dataless/IRISDMC-%s.1P.dataless' % (ista.code))
        sta = newsta.networks[0].stations[0]
    GUA.stations.append(sta)
inv_all.networks.append(GUA)

# AI, which requires other channel names and renaming
inv = client.get_stations(network='AI',station='DSPA',channel='MH?',level='response',\
            minlatitude=-55.5, maxlatitude=-43.1,\
            minlongitude=-76.3, maxlongitude=-65.3)

AI_net = inv.networks[0].select(station='X').copy()
sta = inv.networks[0].select(station='DSPA').stations[0]
ch_list = deepcopy(inv.networks[0].stations[0].channels)
for ch in ch_list:  # the EQ data for this station (MH*) will be resampled to 1Hz to match LH* 
    ch.code = 'L' + ch.code[1:]
    ch.sample_rate = SampleRate(1.)
    ch.location_code = '00'
    sta.channels.append(ch)
AI_net.stations.append(sta)
inv_all.networks.append(AI_net)

# YJ, which requires creating LHE/LHN for some of the stations
yj_rotate = ['IBJ01','IMG01','MEL01','RMB02']
inv = client.get_stations(network='YJ',channel='LH?',level='response',\
            minlatitude=-55.5, maxlatitude=-43.1,\
            minlongitude=-76.3, maxlongitude=-65.3)
YJ_net = inv.networks[0].select(station='X').copy()
for ista in inv.networks[0].stations:
    sta = inv.networks[0].select(station=ista.code).stations[0]
    if ista.code in yj_rotate:
        l2 = inv.select(station=ista.code,channel='LH2').networks[0].stations[0].channels[0].copy()
        l3 = inv.select(station=ista.code,channel='LH3').networks[0].stations[0].channels[0].copy()
        l2.code = 'LHE'; l2.azimuth = Azimuth(90.0)
        l3.code = 'LHN'; l3.azimuth = Azimuth(0.0)
        sta.channels.append(l2)
        sta.channels.append(l3)
    YJ_net.stations.append(sta)
inv_all.networks.append(YJ_net)


# EN
E_net = Network(code='EN',stations=[],\
          description='ENAP from Rodrigo',\
          start_date=UTCDateTime(2019,3,1))
E_file = glob('seed/dataless/ENAP*.xml')
for f in E_file:
    nsta = read_inventory(f)
    for ch in nsta.networks[0].stations[0].channels:
        ch.code = 'L' + ch.code[1:]
    E_net.stations.append(nsta.networks[0].stations[0])
inv_all.networks.append(E_net)

# write everything!
inv_all.write('seed/dataless/combined.xml',format='stationxml')

os.system('java -jar ../other_code/stationxml-seed-converter-2.0.10-SNAPSHOT.jar --input seed/dataless/combined.xml --output seed/dataless/combined.dataless')

