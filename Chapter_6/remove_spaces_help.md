# Remove spaces in gene_presence_absence.csv file
```cut -f1 -d `,'  gene_presence_absence.csv | grep ` ' # finds any spaces ```
# and then:
```sed 's/hdl /hdl_/g' gene_presence_absence.csv > gene_presence_absence_NEW.csv # removes any spaces ```
# Note the 'hdl' is the gene which had spaces in the name
