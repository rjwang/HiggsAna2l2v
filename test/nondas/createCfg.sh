#!bin/sh

logo=$1 #MC8TeV_FermionWIMP_M500_D5_FullSim

cp COPYME_cfg.py mc_${logo}_cfg.py
sed -i "s/COPYME/$logo/g" mc_${logo}_cfg.py

