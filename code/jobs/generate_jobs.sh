#/usr/bin/env bash

# generate the job for latedays
input=$1
threads=$2

if [ ${#} -ne 2 ]; then
  echo "Usage: $0 <input> <threads>"
else
  curdir=`pwd`
  curdir=${curdir%/jobs}
  sed "s:PROGDIR:${curdir}:g" example.job.template > tmp1.job
  sed "s:INPUT:${input}:g" tmp1.job > tmp2.job
  sed "s:THREADS:${threads}:g" tmp2.job > ${USER}_${input}_${threads}.job
  rm -f tmp1.job tmp2.job
fi