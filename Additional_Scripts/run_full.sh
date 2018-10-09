#!/usr/bin/env bash

python Simulation1.py scvi &> ../Sim1/scvi.log
python Simulation1.py others &> ../Sim1/others.log

python Simulation2.py scvi &> ../Sim2/scvi.log
python Simulation2.py others &> ../Sim2/others.log

python Simulation3.py scvi &> ../Sim3/scvi.log
python Simulation3.py others &> ../Sim3/others.log

#python Easycase1.compare.py scvi &> ../Easy1/scvi.log
#python Easycase1.compare.py others &> ../Easy1/others.log

#python Tech1.py scvi &> ../Tech1/scvi.log
#python Tech1.py others &> ../Tech1/others.log

python Easycase2.compare.py scvi &> ../Macosko_Regev/scvi.log
python Easycase2.compare.py others &> ../Macosko_Regev/others.log

python Tech2.py scvi &> ../Zeng/scvi.log
python Tech2.py others &> ../Zeng/others.log

python Simulation1.py scvi &> ../Sim1/scvi.log
python Simulation2.py scvi &> ../Sim2/scvi.log
python Simulation3.py scvi &> ../Sim3/scvi.log
python Easycase1.compare.py scvi &> ../Easy1/scvi.log


