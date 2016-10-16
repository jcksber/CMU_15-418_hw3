#/usr/bin/env bash
# generate jobs in batch

threads=(1 4 16 64 128 240) # The number of threads 
inputs=(easy_4096.txt) # The name of input files

rm -f *.jobs

for f in ${inputs[@]}
do
    for t in ${threads[@]}
    do
        ./generate_jobs.sh $f $t
    done
done
