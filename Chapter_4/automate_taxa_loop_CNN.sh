for i in $(cat AB_taxa.txt); 
  do python3 Taxa_CNN_loop.py -test_file "taxa_CNN_file/${i}_eggnog_taxa_test_file.csv" -training_data "taxa_CNN_file/${i}_not_eggnog_taxa_train_file.csv" -result_file "${i}"; 
done
