#!/bin/bash

/scratch/grudorff/opt/usr/bin/parallel -j 10 /scratch/grudorff/rasmus/python/msc-rasmus/cp2k-forces/analysis.sh -i /scratch/grudorff/hfxforcetest/run-0/hematite_localLog_p{}.out ::: {1..767}
