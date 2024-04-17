# Eggnog_output_to_arff
Translates eggnog output to .arff file to be used in weka, specifically the most general 'Cog' or gene family (absence or number of copies of each gene). 
I merge with some other data to see the relationship between gene families and phenotype of my strains.

This is currently very specific to my work but menu can be altered for your needs. (Still working on writing code to be more universal).



usage: eggnog_output_to_arff.py [-h] -file FILE -o OUTPUT_FILE -c_attr_name
                                CAT_ATTR_NAME -catagory1 CATAGORY1 -catagory2
                                CATAGORY2 [-catagory3 CATAGORY3]

Run eggnog to arff tool parameters

optional arguments:
  -h, --help            show this help message and exit
  -file FILE            eggnog output file
  -o OUTPUT_FILE        Name of output file ext .arff
  -c_attr_name CAT_ATTR_NAME
                        Catagorical variable column title
  -catagory1 CATAGORY1  Catagory 1 e.g. Susceptible
  -catagory2 CATAGORY2  Catagory 2 e.g. Resistant
  -catagory3 CATAGORY3  Catagory 3 e.g. Intermediate


  Note: If you are using the eggnog file on the original model (dot files apply decision tree) it needs to be GeneFamily- not Cog.
