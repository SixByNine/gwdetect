#!/bin/bash


set -- `meanvar.py onstats 4`
mean=$1
set -- `meanvar.py offstats 4`
sig=$3

sn=`echo $mean $sig | awk '{print $1/$2}'`

echo "$sn (4)"

set -- `meanvar.py onstats 1`
mean=$1
set -- `meanvar.py offstats 1`
sig=$3

sn=`echo $mean $sig | awk '{print $1/$2}'`

echo "$sn (1)"

