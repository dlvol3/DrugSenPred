#!/bin/bash -l
#OAR -n multSNE
#OAR -l nodes=1/core=20,walltime=52:00:00

module use /opt/easybuild/modules/all
module load R/3.4.0-foss-2017a-X11-20170314

cd /mnt/pcpnfs/homedirs/yzhang/CCLEGDSC_RF/scripts/

R CMD BATCH GDSCtrain_CCLEtest_remote.r
