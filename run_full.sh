for f in Easycase1.compare.py Simulation2.py Simulation1.py Simulation3.py Tech1.py
do 
for method in vae scanvi scanvi1 scanvi2 Seurat Combat
do
python $f $method 
done >> $f.res.txt
done

for f in Easycase2.compare.py Tech2.py
do
for method in vae scanvi scanvi1 scanvi2 readSeurat Combat
do
python $f $method 
done >> $f.res.txt
done

for method in vae scanvi scanvi1 scanvi2 Seurat Combat
do
python Easycase1.compare.py $method  small
done >> $f.small.res.txt

for method in vae scanvi scanvi1 scanvi2 readSeurat Combat
do
python Easycase1.compare.py $method  large
done >> $f.large.res.txt

