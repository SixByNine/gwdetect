#!/bin/bash
amp=1e-14
nthread=8
. ./settings

rm real.*/hd.plot real.*/ff.plot fpy*

for trial in `seq 0 $(($nreal-1))` ; do
   if [[ -e real.$trial/result.dat ]] ; then
	  if [[ ! -e pres.$trial ]] ;then

		 while [[ `jobs -r | wc -l` -ge $nthread ]] ; do
			n=`jobs -r | wc -l`
			sleep 0.1
		 done

		 echo -ne "\rReal $trial  "
		 detect_GWB.py real.$trial/GW.sum real.$trial $amp -Lfcvm -ffactors $* > fpy.$trial &
	  fi
   fi
done

echo "Wait for threads to end"
wait
