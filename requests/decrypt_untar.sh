#!/bin/bash
#

## decrypt data with openssl, untar, and get rid of encrypted files and archives
# give decryption password as a command line argument

if [ "$1" == "" ]; then
   echo 'please supply a password for decryption'
   exit 0
fi

#seeddir=../seed/1P
seeddir=../seed/second_half_1P


cd $seeddir

for encr in  *'.openssl'; do
  /usr/bin/openssl enc -d -des-cbc -salt -in ${encr} -out ${encr%.openssl} -pass pass:$1
  tar -xf ${encr%.openssl}
  rm -f ${encr}
  rm -f ${encr%.openssl}

done
