#/usr/bin/env bash
# generate jobs in batch

threads=(1 2 4) # The number of threads 
inputs=(easy_16384.txt) # The name of input files

rm -f *.jobs

for f in ${inputs[@]}
do
    for t in ${threads[@]}
    do
        ./generate_jobs.sh $f $t
    done
done
