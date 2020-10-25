import numpy as np
from obspy import read_inventory, UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core.inventory.util import Azimuth, SampleRate
from glob import glob
from copy import deepcopy
import os, sys

####
# read stationxml files for instruments, get lists of the ones we actually use,
# combine into one stationxml, write, convert to dataless seed for rdseed, cross fingers
####

inv_all = Inventory(networks=[],source='IRIS_ENAP')

# read C, find insts, add just those
C_read = read_inventory('seed/dataless/C.xml')
C_file = glob('seed/dataless/IRISDMC*.C.dataless')
C_inst = [e.split('/')[-1].split('-')[-1].split('.')[0] for e in C_file]

C_net = C_read.networks[0].select(station='X').copy()  # just the network
for i in C_inst:
    sta = C_read.networks[0].select(station=i).stations[0]
    C_net.stations.append(sta)
inv_all.networks.append(C_net)

del C_read  # cleanup

# read C1, find insts, add just those
C1_read = read_inventory('seed/dataless/C1.xml')
C1_file = glob('seed/dataless/IRISDMC*.C1.dataless')
C1_inst = [e.split('/')[-1].split('-')[-1].split('.')[0] for e in C1_file]

if 'MG01' not in C1_inst:
    C1_inst.append('MG01')

C1_net = C1_read.networks[0].select(station='X').copy()  # just the network
for i in C1_inst:
    sta = C1_read.networks[0].select(station=i).stations[0]
    C1_net.stations.append(sta)
inv_all.networks.append(C1_net)

del C1_read  # cleanup

# read G, find insts, add just those
G_read = read_inventory('seed/dataless/G.xml')
G_file = glob('seed/dataless/IRISDMC*.G.dataless')
G_inst = [e.split('/')[-1].split('-')[-1].split('.')[0] for e in G_file]

G_net = G_read.networks[0].select(station='X').copy()  # just the network
for i in G_inst:
    sta = G_read.networks[0].select(station=i).stations[0]
    G_net.stations.append(sta)
inv_all.networks.append(G_net)

del G_read  # cleanup

# read AI, find the one station, add it
AI_read = read_inventory('seed/dataless/AI.dataless')

AI_net = AI_read.networks[0].select(station='X').copy()  # just the network
sta = AI_read.networks[0].select(station='DSPA').stations[0]
ch_list = deepcopy(sta.channels[-3:])  # we happen to know these are MH*
for ch in ch_list:
    ch.code = 'L' + ch.code[1:]
    ch.sample_rate = SampleRate(1.)  # for resampled LH* data, when it gets rewritten
    ch.location_code = '00'  # for simplicity
    sta.channels.append(ch)

AI_net.stations.append(sta)
inv_all.networks.append(AI_net)

del AI_read  # cleanup

# read XB, find insts, add just those
XB_read = read_inventory('seed/dataless/XB.1997-1999.dataless')
XB_file = glob('seed/dataless/IRISDMC*.XB.dataless')
XB_inst = [e.split('/')[-1].split('-')[-1].split('.')[0] for e in XB_file]

XB_net = XB_read.networks[0].select(station='X').copy()  # just the network
for i in XB_inst:
    sta = XB_read.networks[0].select(station=i).stations[0]
    XB_net.stations.append(sta)
inv_all.networks.append(XB_net)

del XB_read  # cleanup

# read YN, find insts, add just those  # NONE OF THESE ARE USABLE so skip the network
#YN_read = read_inventory('seed/dataless/YN.1999-2004.dataless')
#YN_file = glob('seed/dataless/IRISDMC*.YN.dataless')
#YN_inst = [e.split('/')[-1].split('-')[-1].split('.')[0] for e in YN_file]

#YN_net = YN_read.networks[0].select(station='X').copy()  # just the network
#for i in YN_inst:
#    sta = YN_read.networks[0].select(station=i).stations[0]
#    YN_net.stations.append(sta)
#inv_all.networks.append(YN_net)

#del YN_read  # cleanup

# read 1P, find insts, add just those
P_read = read_inventory('seed/dataless/1P.xml')
P_file = glob('seed/dataless/IRISDMC*.1P.dataless')
P_inst = [e.split('/')[-1].split('-')[-1].split('.')[0] for e in P_file]

P_net = P_read.networks[0].select(station='X').copy()  # just the network
bad_sta = []
for i in P_inst:
    try:
        sta = P_read.networks[0].select(station=i).stations[0]
    except IndexError:
        bad_sta.append(i)
        continue
    P_net.stations.append(sta)

# add the stations that maybe changed names?
for i in bad_sta:
    if i == 'BRMY':  # rename the renamed stations
        sta2 = 'BRNY'
        sta = P_read.networks[0].select(station=sta2).stations[0]
        sta.code = 'BRMY'
        P_net.stations.append(sta)
    elif i == 'PTTY':
        sta2 = 'CH14'
        sta = P_read.networks[0].select(station=sta2).stations[0]
        sta.code = 'PTTY'
        P_net.stations.append(sta)
    else:
        os.system('java -jar ../other_code/stationxml-seed-converter-2.0.10-SNAPSHOT.jar --input seed/dataless/IRISDMC-%s.1P.dataless --output seed/dataless/%s.1P.xml' % (i,i))
        temp_inv = read_inventory('seed/dataless/%s.1P.xml' % i)
        P_net.stations.append(temp_inv.networks[0].stations[0])

inv_all.networks.append(P_net)

del P_read  # cleanup

# read YJ, find insts, add just those
yj_rotate = ['IBJ01','IMG01','MEL01']  # LHZ/2/3 stations; add extra 'rotated' channels
Y_read = read_inventory('seed/dataless/YJ.xml')
Y_file = glob('seed/dataless/IRISDMC*.YJ.dataless')
Y_inst = [e.split('/')[-1].split('-')[-1].split('.')[0] for e in Y_file]

Y_net = Y_read.networks[0].select(station='X').copy()  # just the network
bad_sta = []
for i in Y_inst:
    sta = Y_read.networks[0].select(station=i).stations[0]
    if i in yj_rotate:
        # copy LH2 and LH3 channels
        l2 = Y_read.select(station=i,channel='LH2').networks[0].stations[0].channels[0].copy()
        l3 = Y_read.select(station=i,channel='LH3').networks[0].stations[0].channels[0].copy()
        l2.code = 'LHE'; l2.azimuth = Azimuth(90.0)
        l3.code = 'LHN'; l3.azimuth = Azimuth(0.0)
        sta.channels.append(l2)
        sta.channels.append(l3)
    Y_net.stations.append(sta)
inv_all.networks.append(Y_net)

del Y_read  # cleanup


# deal with ENAP stations
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
