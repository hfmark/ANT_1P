import numpy as np
from datetime import datetime, timedelta
import os, sys


# XB and YN (SEPA stations)

stalst = '/P/hmark/ANT_1P/verysouth_sta.lst'  # network  station
net,sta = np.loadtxt(stalst,dtype=(str,str),unpack=True)

request_name = 'Hannah Mark'
request_email = 'hmark@wustl.edu'

day1 = datetime(1997,1,6)
dayN = datetime(2001,12,31)

dayR = day1
while dayR < dayN:
	BT = dayR

	label = 'SP_LH_%s' % (dayR.strftime('%Y.%b.%d'))

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
			mailinfo.write('1 LH?\n')
	mailinfo.close()
	print('Start sending mail')
	os.system("cat %s | mail -s 'Requesting Data' breq_fast@iris.washington.edu" % mainname)
	print('Mail sent')

