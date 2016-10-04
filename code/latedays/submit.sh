#/usr/bin/env bash
rm -f $(whoami)_*
for job in ../jobs/*.job
do
    qsub $job
done
