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

# deal with AY03: 120P up until 2 Nov 2017 (or 11 Feb?), T40 after (T40 response is on DMC)
i2 = _read_seed('seed/dataless/IRISDMC-AY03.C1.dataless')
sta = inv.select(station='AY03').networks[0].stations[0]
ch_list_add = deepcopy(i2.networks[0].stations[0].channels)
ch_list_init = deepcopy(sta.channels)

for ch in ch_list_init:
    ch.start_date = UTCDateTime('2017-11-02')  # or is it 2017-02-11???

for ch in ch_list_add:
    ch.end_date = UTCDateTime('2017-11-02')
    ch_list_init.append(ch)

# somehow find AY03 in the overall inv and replace the channel list
codes = np.array([e.code for e in inv.networks])
inet = np.where(codes == 'C1')[0][0]
nms = np.array([e.code for e in inv.networks[inet].stations])
ista = np.where(nms == 'AY03')[0][0]
inv.networks[inet].stations[ista].channels = ch_list_init

# add all these networks to the main inventory
for net in inv.networks:
    inv_all.networks.append(net)


# 1P, which is read from Patrick's file (at least until DMC is updated with horizons)
inv = _read_seed('seed/dataless/1P.20.patadb.202104092300.dataless')
# will have junk of lots of channels not being used but oh well
inv_all.networks.append(inv.networks[0])

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

esta = read_inventory('seed/dataless/ENAP-ANMA.EN.xml')
for ch in esta.networks[0].stations[0].channels:
    ch.code = 'L' + ch.code[1:]
inv_all.networks[-1].stations.append(esta.networks[0].stations[0])

# write everything!
inv_all.write('seed/dataless/combined.xml',format='stationxml')

os.system('java -jar ../other_code/stationxml-seed-converter-2.0.10-SNAPSHOT.jar --input seed/dataless/combined.xml --output seed/dataless/combined.dataless')

