# Pipe line commands for converting NCBI AMRFinder output to gene counts:


#### Its a good idea to make a unique list of all the genes to help you identify gene hits:
#### Gene hits are in the 9th column of RGI output txt files, you can use the following command:

``` 
cut -f6 *.fna.txt |tr -d "'" | sed '1d' | sort | uniq > NCBI_genes.txt 

```
#### Cut the gene hit columns for each genome:

```
for i in *.fasta.txt; do cut -f6 $i | tr -d "'" | sed '1d' > $i.txt; done
```
 #### Make gene count files:
 #### This shell script makes an output file of the genes and the gene count e.g. if it has a 1 hit, it would be "gene_name 1" 

``` 
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --partition=medpri
#SBATCH --job-name=gene_count_RGI
#SBATCH --array=1-167%100
echo $SLURM_ARRAY_TASK_ID
echo $SLURM_JOB_ID
cd /mnt/scratch2/users/40309916/RGI_output


for i in ls *.txt.txt; 
do 
sh NCBI_count_script.sh $i > $i.NCBIcount.txt;
done
```
#### NCBI_count_script.sh : (to be used in script above)
```
cat $1 | awk 'BEGIN{ 
    
    while(( getline line < "NCBI_genes.txt") > 0 ) { 
           array[line]=0; 
           } 

    };{
    
    if( $1 in array)
        {
        array[$1] = array[$1]+1;
        }
    else
        print "ERROR $1 not in array";
    
    }END{ 
    
    for(gene in array)
        {
        print gene"\t"array[gene]
        }
    
}'
```

#### This next script combines all data so you can analyse, data as one dataframe: catfiles.sh
```
#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --partition=hipri
#SBATCH --job-name=cat_files
cd /mnt/scratch2/users/40309916/NCBI_AMRFinder
for i in *count.txt;
do
name=`echo $i | cut -d'_' -f2`;
sed "s/^/$name /g" < $i;
done > myinputfile.txt
```
