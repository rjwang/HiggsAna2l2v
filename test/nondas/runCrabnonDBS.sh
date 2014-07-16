#!/bin/sh


echo "+++++++++++++++++++++++++++++++++++"
echo "          runCrabnonDBS.sh         "
echo "            Renjie Wang            "
echo "        renjie.wang@cern.ch        "
echo "+++++++++++++++++++++++++++++++++++"

i="$1"

if [ "$i" == "help" ]; then
echo ">>>usage: testCrabnonDBS.sh <job index> [<max events>]"
  exit 0;
fi

if [ "$i" == "" ]; then
echo ">>>Error: missing job index"
  exit 1;
fi

echo ">>> job number: #$i "
python generateSingleFileCfg.py -N $i
cmsRun -j $RUNTIME_AREA/crab_fjr_$i.xml -p pset.py


echo "+++++++++++++++++++++++++++++++++++"
echo "          runCrabnonDBS.sh         "
echo "+++++++++++++++++++++++++++++++++++"

