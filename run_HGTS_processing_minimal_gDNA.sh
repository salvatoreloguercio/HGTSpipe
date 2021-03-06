
jobname=`echo R_$1 | cut -c1-15`

echo "#PBS -N ${jobname}"
echo '#PBS -l nodes=1:ppn=8'
echo '#PBS -l mem=30gb'
echo '#PBS -l walltime=240:00:00'
echo '#PBS -l cput=9600:00:00'
echo '#PBS -j oe'
echo '#PBS -m ea'
echo '#PBS -M ekleiman@scripps.edu'
echo
echo 'cd $PBS_O_WORKDIR'
echo 'module load R'
echo
echo "Rscript --no-restore --no-save HGTS_processing_minimal.R -i $1 -e $2 -c $3 -t gDNA -d $4 -l $5 -u $6 -n $7"
echo
