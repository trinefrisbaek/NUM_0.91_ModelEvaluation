#!/bin/sh
i=$1
echo $'nohup matlab -nodisplay -r "baserunWatercolumn('$1')">&output'$1'.txt&' >> myjob.pbs
