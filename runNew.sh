#!/bin/bash
make main.exe -j
if (($?)); then
exit -1
fi

bash backupSource.sh
export OMP_NUM_THREADS=4
./main.exe