import numpy as np
from datetime import datetime, timedelta
import os, sys

# permanent stations, prior to 1P

#stalst = '/P/hmark/ANT_1P/requests/permanent_sta.lst'  # network  station
#stalst = '/P/hmark/ANT_1P/requests/extra_permanent.lst'  # network  station
stalst = '/P/hmark/ANT_1P/requests/AI_special.lst'  # network  station
net,sta = np.loadtxt(stalst,dtype=(str,str),unpack=True)

request_name = 'Hannah Mark'
request_email = 'hmark@wustl.edu'

#day1 = datetime(2010,7,16)  # all dates for at least 2 permanent stations up to start of 1P
#dayN = datetime(2018,11,4)
day1 = datetime(2010,7,16)  # all dates for at least 2 permanent stations INCLUDING 1P
dayN = datetime(2019,11,1) 

dayR = day1
while dayR < dayN:
	BT = dayR

	#label = '1P_LH_%s' % (dayR.strftime('%Y.%b.%d'))
	label = '1P_MH_%s' % (dayR.strftime('%Y.%b.%d'))

	dayR = dayR + timedelta(days=1)
	ET = dayR

	sss = '''.NAME %s
.INST Washington University in St. Louis
.MAIL 1 Brookings Drive
.EMAIL %s
.PHONE
.FAX
.MEDIA: Electronic (FTP)
.LABEL %s
.END
''' % (request_name, request_email, label)

	mainname = 'mail.tem'
	mailinfo = open(mainname,'w')
	mailinfo.write(sss)
	with open(stalst,'r') as stal:
		for sline in stal.readlines():
			slineinfo = sline.split()
			mailinfo.write('%s %s %s %s %s %s %s %s ' % \
				(slineinfo[1], slineinfo[0], BT.year, BT.month, \
				BT.day, BT.hour, BT.minute, BT.second))
			mailinfo.write('%s %s %s %s %s %s ' % (ET.year, ET.month, ET.day, \
								    ET.hour, ET.minute, ET.second))
			#mailinfo.write('1 LH?\n')
			mailinfo.write('1 MH?\n')
	mailinfo.close()
	print('Start sending mail')
	os.system("cat %s | mail -s 'Requesting Data' breq_fast@iris.washington.edu" % mainname)
	print('Mail sent')

