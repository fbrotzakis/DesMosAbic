for ((c=0;c<61;c++));do
echo $c
dir=r"$c"
echo $dir
cd $dir

cp ../scripts/make_pdb_files.py .
cp ../scripts/design.py .
cp ../scripts/submit.job .
rm *pdb 
rm -r output
sbatch submit.job 
#python make_pdb_files.py
cd -
done

