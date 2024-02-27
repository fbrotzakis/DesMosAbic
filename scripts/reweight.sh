rm JUMPTIMES_MT JUMPTIMES_WT weights.dat 
for ((c=0;c<61;c++));do
echo $c
dir=r"$c"
echo $dir

cd $dir/output/
cp ../../scripts/plumed_analysis.dat .
cp ../../scripts/structure.pdb .
cp ../../scripts/protein*.pdb .
plumed driver --mf_xtc ../run.xtc --plumed plumed_analysis.dat
cp ../../scripts/reweight.py .
python reweight.py
cd -
done

python cal_cdf.py
sort -g JUMPTIMES_WT |awk '{print sum+=$1}' > CDF_WT.dat
sort -g JUMPTIMES_MT |awk '{print sum+=$1}' > CDF_MT.dat
python cdf_plt.py
python DKL.py
