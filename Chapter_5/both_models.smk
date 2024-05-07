singularity: "./ldillon_amrwebsite.sif"

GENOMES, = glob_wildcards("static/files/{genome}.fna")

rule all:
    input:
        expand("static/files/results/{genome}.txt", genome=GENOMES)


rule prodigal:
    input:
        "static/files/{genome}.fna"
    output:
        temp("static/files/prodigal_output/{genome}.faa")
    singularity:
        "ldillon_amrwebsite.sif"
    shell:
        '''
        prodigal -i {input} -a {output} -f gff
        '''

rule diamond:
    input:
        faa="static/files/prodigal_output/{genome}.faa",
        db="data/eggnog_proteins.dmnd"
    output:
        temp("static/files/diamond_annotation/{genome}_Diamond.faa")
    shell:
        '''
        diamond blastp -d {input.db} -q {input.faa} --more-sensitive --threads 8 -e 0.001000 -o {output} --top 3 --query-cover 0 --subject-cover 0
        '''

rule Edit_diamond:
    input:
        "static/files/diamond_annotation/{genome}_Diamond.faa"
    output:
        temp("static/files/diamond_annotation/{genome}_DiamondReduced.faa")
    run:
        shell("""cat {input} | awk 'BEGIN{{prev = ""}};{{if ($1 != prev) {{ prev=$1; print $1"\t"$2"\t"$11"\t"$12 }} }}' > {output}""")


rule eggnog:
    input:
        faa="static/files/diamond_annotation/{genome}_DiamondReduced.faa",
        db="data/"
    output:
        temp("static/files/{genome}.emapper.annotations")
    params:
        prefix="{genome}"
    shell:
        """
        emapper.py --data_dir {input.db} --annotate_hits_table {input.faa}  --no_file_comments --override --output {params.prefix} --cpu 8;
        mv {params.prefix}.emapper.annotations static/files/
        """

rule eggnog_2_arff:
    input:
        "static/files/{genome}.emapper.annotations"
    output:
        temp("static/files/decision_tree/{genome}.arff")
    shell:
        """
        python eggnog_output_to_arff.py -file {input} -o {output} -c_attr_name phenotype -catagory1 Susceptible -catagory2 Resistant
        """

rule genome_names:
    input:
        "static/files/{genome}.fna"
    output:
        temp("static/files/decision_tree/{genome}_names.txt")
    run:
        shell("""ls {input} > {output}; echo >> {output}""")


rule run_c_program: 
    input:
        dot_ami="dot_files/amikacin_EGGNOG_NEW.arff.dot",
        dot_amo="dot_files/amoxicillin_EGGNOG_NEW.arff.dot",
        dot_amp="dot_files/ampicillin_EGGNOG_NEW.arff.dot",
        dot_az="dot_files/aztreonam_EGGNOG_NEW.arff.dot",
        dot_cefe="dot_files/cefepime_EGGNOG_NEW.arff.dot",
        dot_ceft="dot_files/ceftriaxone_EGGNOG_NEW.arff.dot",
        dot_chlo="dot_files/chloramphenicol_EGGNOG_NEW.arff.dot",
        dot_cip="dot_files/ciprofloxacin_EGGNOG_NEW.arff.dot",
        dot_clin="dot_files/clindamycin_EGGNOG_NEW.arff.dot",
        dot_col="dot_files/colistin_EGGNOG_NEW.arff.dot",
        dot_dor="dot_files/doripenem_EGGNOG_NEW.arff.dot",
        dot_ert="dot_files/ertapenem_EGGNOG_NEW.arff.dot",
        dot_eryt="dot_files/erythromycin_EGGNOG_NEW.arff.dot",
        dot_fos="dot_files/fosfomycin_EGGNOG_NEW.arff.dot",
        dot_gen="dot_files/gentamicin_EGGNOG_NEW.arff.dot",
        dot_imi="dot_files/imipenem_EGGNOG_NEW.arff.dot",
        dot_lev="dot_files/levofloxacin_EGGNOG_NEW.arff.dot",
        dot_mer="dot_files/meropenem_EGGNOG_NEW.arff.dot",
        dot_mox="dot_files/moxifloxacin_EGGNOG_NEW.arff.dot",
        dot_nit="dot_files/nitrofurantoin_EGGNOG_NEW.arff.dot",
        dot_tet="dot_files/tetracycline_EGGNOG_NEW.arff.dot",
        dot_tig="dot_files/tigecycline_EGGNOG_NEW.arff.dot",
        dot_tob="dot_files/tobramycin_EGGNOG_NEW.arff.dot",
        genome_names="static/files/decision_tree/{genome}_names.txt",
        arff_file="static/files/decision_tree/{genome}.arff",
    output:
        paths=temp("static/files/decision_tree/{genome}_names.txt.paths.txt"),
        predictions=temp("static/files/decision_tree/{genome}_names.txt.predicitions.txt"),
        predictions_ami=temp("static/files/decision_tree/{genome}_ami.txt.predicitions.txt"),
        predictions_amo=temp("static/files/decision_tree/{genome}_amo.txt.predicitions.txt"),
        predictions_amp=temp("static/files/decision_tree/{genome}_amp.txt.predicitions.txt"),
        predictions_az=temp("static/files/decision_tree/{genome}_az.txt.predicitions.txt"),
        predictions_cefe=temp("static/files/decision_tree/{genome}_cefe.txt.predicitions.txt"),
        predictions_ceft=temp("static/files/decision_tree/{genome}_ceft.txt.predicitions.txt"),
        predictions_chlo=temp("static/files/decision_tree/{genome}_chlo.txt.predicitions.txt"),
        predictions_cip=temp("static/files/decision_tree/{genome}_cip.txt.predicitions.txt"),
        predictions_clin=temp("static/files//decision_tree/{genome}_clin.txt.predicitions.txt"),
        predictions_col=temp("static/files/decision_tree/{genome}_col.txt.predicitions.txt"),
        predictions_dor=temp("static/files/decision_tree/{genome}_dor.txt.predicitions.txt"),
        predictions_ert=temp("static/files/decision_tree/{genome}_ert.txt.predicitions.txt"),
        predictions_eryt=temp("static/files/decision_tree/{genome}_eryt.txt.predicitions.txt"),
        predictions_fos=temp("static/files/decision_tree/{genome}_fos.txt.predicitions.txt"),
        predictions_gen=temp("static/files/decision_tree/{genome}_gen.txt.predicitions.txt"),
        predictions_imi=temp("static/files/decision_tree/{genome}_imi.txt.predicitions.txt"),
        predictions_lev=temp("static/files/decision_tree/{genome}_lev.txt.predicitions.txt"),
        predictions_mer=temp("static/files/decision_tree/{genome}_mer.txt.predicitions.txt"),
        predictions_mox=temp("static/files/decision_tree/{genome}_mox.txt.predicitions.txt"),
        predictions_nit=temp("static/files/decision_tree/{genome}_nit.txt.predicitions.txt"),
        predictions_tet=temp("static/files/decision_tree/{genome}_tet.txt.predicitions.txt"),
        predictions_tig=temp("static/files/decision_tree/{genome}_tig.txt.predicitions.txt"),
        predictions_tob=temp("static/files/decision_tree/{genome}_tob.txt.predicitions.txt"),
    shell:
        """
        my_apply_decision_tree {input.dot_ami} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_ami};
        my_apply_decision_tree {input.dot_amo} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_amo};
        my_apply_decision_tree {input.dot_amp} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_amp};
        my_apply_decision_tree {input.dot_az} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_az};
        my_apply_decision_tree {input.dot_cefe} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_cefe};
        my_apply_decision_tree {input.dot_ceft} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_ceft};
        my_apply_decision_tree {input.dot_chlo} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_chlo};
        my_apply_decision_tree {input.dot_cip} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_cip};
        my_apply_decision_tree {input.dot_clin} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_clin};
        my_apply_decision_tree {input.dot_col} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_col};
        my_apply_decision_tree {input.dot_dor} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_dor};
        my_apply_decision_tree {input.dot_ert} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_ert};
        my_apply_decision_tree {input.dot_eryt} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_eryt};
        my_apply_decision_tree {input.dot_fos} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_fos};
        my_apply_decision_tree {input.dot_gen} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_gen};
        my_apply_decision_tree {input.dot_imi} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_imi};
        my_apply_decision_tree {input.dot_lev} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_lev};
        my_apply_decision_tree {input.dot_mer} {input.arff_file} {input.genome_names};
        tail -1  {output.predictions} > {output.predictions_mer};
        my_apply_decision_tree {input.dot_mox} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_mox};
        my_apply_decision_tree {input.dot_nit} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_nit};
        my_apply_decision_tree {input.dot_tet} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_tet};
        my_apply_decision_tree {input.dot_tig} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_tig};
        my_apply_decision_tree {input.dot_tob} {input.arff_file} {input.genome_names};
        tail -1 {output.predictions} > {output.predictions_tob};
        """

rule combine_files:
    input:
        predictions_cip="static/files/decision_tree/{genome}_cip.txt.predicitions.txt",
        predictions_ami="static/files/decision_tree/{genome}_ami.txt.predicitions.txt",
        predictions_amo="static/files/decision_tree/{genome}_amo.txt.predicitions.txt",
        predictions_amp="static/files/decision_tree/{genome}_amp.txt.predicitions.txt",
        predictions_az="static/files/decision_tree/{genome}_az.txt.predicitions.txt",
        predictions_cefe="static/files/decision_tree/{genome}_cefe.txt.predicitions.txt",
        predictions_ceft="static/files/decision_tree/{genome}_ceft.txt.predicitions.txt",
        predictions_chlo="static/files/decision_tree/{genome}_chlo.txt.predicitions.txt",
        predictions_cip="static/files/decision_tree/{genome}_cip.txt.predicitions.txt",
        predictions_clin="static/files/decision_tree/{genome}_clin.txt.predicitions.txt",
        predictions_col="static/files/decision_tree/{genome}_col.txt.predicitions.txt",
        predictions_dor="static/files/decision_tree/{genome}_dor.txt.predicitions.txt",
        predictions_ert="static/files/decision_tree/{genome}_ert.txt.predicitions.txt",
        predictions_eryt="static/files/decision_tree/{genome}_eryt.txt.predicitions.txt",
        predictions_fos="static/files/decision_tree/{genome}_fos.txt.predicitions.txt",
        predictions_gen="static/files/decision_tree/{genome}_gen.txt.predicitions.txt",
        predictions_imi="static/files/decision_tree/{genome}_imi.txt.predicitions.txt",
        predictions_lev="static/files/decision_tree/{genome}_lev.txt.predicitions.txt",
        predictions_mer="static/files/decision_tree/{genome}_mer.txt.predicitions.txt",
        predictions_mox="static/files/decision_tree/{genome}_mox.txt.predicitions.txt",
        predictions_nit="static/files/decision_tree/{genome}_nit.txt.predicitions.txt",
        predictions_tet="static/files/decision_tree/{genome}_tet.txt.predicitions.txt",
        predictions_tig="static/files/decision_tree/{genome}_tig.txt.predicitions.txt",
        predictions_tob="static/files/decision_tree/{genome}_tob.txt.predicitions.txt"
    output:
        final=temp("static/files/decision_tree/{genome}_FINAL.txt.DT.predicitions.txt")
    run:
        shell("""echo "Antibiotic	Genome_name	Genome_number	Node_number	Node_label" > {output.final}; awk '{{print "Amikacin", $0}}' {input.predictions_ami} >> {output.final}; awk '{{print "Amoxicillin", $0}}' {input.predictions_amo} >> {output.final}; awk '{{print "Ampicillin", $0}}' {input.predictions_amp} >> {output.final}; awk '{{print "Aztreonam", $0}}' {input.predictions_az} >> {output.final}; awk '{{print "Cefepime", $0}}' {input.predictions_cefe} >> {output.final}; awk '{{print "Ceftriaxone", $0}}' {input.predictions_ceft} >> {output.final}; awk '{{print "Chloramphenicol", $0}}' {input.predictions_chlo} >> {output.final}; awk '{{print "Ciprofloxacin", $0}}' {input.predictions_cip} >> {output.final}; awk '{{print "Clindamycin", $0}}' {input.predictions_clin} >> {output.final}; awk '{{print "Colistin", $0}}' {input.predictions_col} >> {output.final}; awk '{{print "Doripenem", $0}}' {input.predictions_dor} >> {output.final}; awk '{{print "Ertapenem", $0}}' {input.predictions_ert} >> {output.final}; awk '{{print "Erythromcyin", $0}}' {input.predictions_eryt} >> {output.final}; awk '{{print "Fosfomycin", $0}}' {input.predictions_fos} >> {output.final}; awk '{{print "Gentamicin", $0}}' {input.predictions_gen} >> {output.final}; awk '{{print "Imipenem, $0}}' {input.predictions_imi} >> {output.final}; awk '{{print "Levofloxacin", $0}}' {input.predictions_lev} >> {output.final}; awk '{{print "Meropenem", $0}}' {input.predictions_mer} >> {output.final}; awk '{{print "Moxifloxacin", $0}}' {input.predictions_mox} >> {output.final}; awk '{{print "Nitrofurantoin", $0}}' {input.predictions_nit} >> {output.final}; awk '{{print "Tetracycline", $0}}' {input.predictions_tet} >> {output.final}; awk '{{print "Tigecycline", $0}}' {input.predictions_tig} >> {output.final}; awk '{{print "Tobramycin", $0}}' {input.predictions_tob} >> {output.final}""")

# Now convert the eggNOG output to an appropriate format for python
rule eggnog2cnn:
    input:
         eggnog_output="{genome}.emapper.annotations"
         cnn_ami="Training_models/amikacin_egg_model.h5"
         cnn_amo="Training_models/amoxicillin_egg_model.h5"
         cnn_amp="Training_models/ampicillin_egg_model.h5"
         cnn_az="Training_models/aztreonam_egg_model.h5"
         cnn_cefe="Training_models/cefepime_egg_model.h5"
         cnn_ceft="Training_models/ceftriaxone_egg_model.h5"
         cnn_chlo="Training_models/chloramphenicol_egg_model.h5"
         cnn_cip="Training_models/ciprofloxacin_egg_model.h5"
         cnn_clin="Training_models/clindamycin_egg_model.h5"
         cnn_col="Training_models/colistin_egg_model.h5"
         cnn_dor="Training_models/doripenem_egg_model.h5"
         cnn_ert="Training_models/ertapenem_egg_model.h5"
         cnn_eryt="Training_models/erythromycin_egg_model.h5"
         cnn_fos="Training_models/fosfomycin_egg_model.h5"
         cnn_gen="Training_models/gentamicin_egg_model.h5"
         cnn_imi="Training_models/imipenem_egg_model.h5"
         cnn_lev="Training_models/levofloxacin_egg_model.h5"
         cnn_mer="Training_models/meropenem_egg_model.h5"
         cnn_mox="Training_models/moxifloxacin_egg_model.h5"
         cnn_nit="Training_models/nitrofurantoin_egg_model.h5"
         cnn_tet="Training_models/tetracycline_egg_model.h5"
         cnn_tig="Training_models/tigecycline_egg_model.h5"
         cnn_tob="Training_models/tobramycin_egg_model.h5"
    output:
        output_ami="CNN_predictions_ami.txt"
        output_amo="CNN_predictions_amo.txt"
        output_amp="CNN_predictions_amp.txt"
        output_az="CNN_predictions_az.txt"
        output_cefe="CNN_predictions_cefe.txt"
        output_ceft="CNN_predictions_ceft.txt"
        output_chlo="CNN_predictions_chlo.txt"
        output_cip="CNN_predictions_cip.txt"
        output_clin="CNN_predictions_clin.txt"
        output_col="CNN_predictions_col.txt"
        output_dor="CNN_predictions_dor.txt"
        output_ert="CNN_predictions_ert.txt"
        output_eryt="CNN_predictions_eryt.txt"
        output_fos="CNN_predictions_fos.txt"
        output_gen="CNN_predictions_gen.txt"
        output_imi="CNN_predictions_imi.txt"
        output_lev="CNN_predictions_lev.txt"
        output_mer="CNN_predictions_mer.txt"
        output_mox="CNN_predictions_mox.txt"
        output_nit="CNN_predictions_nit.txt"
        output_tet="CNN_predictions_tet.txt"
        output_tig="CNN_predictions_tig.txt"
        output_tob="CNN_predictions_tob.txt"
    shell:
        """
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_ami} -o {output_ami};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_amo} -o {output_amo};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_amp} -o {output_amp};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_az} -o {output_az};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_cefe} -o {output_cefe};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_ceft} -o {output_ceft};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_chlo} -o {output_chlo};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_cip} -o {output_cip};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_clin} -o {output_clin};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_col} -o {output_col};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_dor} -o {output_dor};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_ert} -o {output_ert};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_eryt} -o {output_eryt};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_fos} -o {output_fos};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_gen} -o {output_gen};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_imi} -o {output_imi};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_lev} -o {output_lev};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_mer} -o {output_mer};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_mox} -o {output_mox};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_nit} -o {output_nit};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_tet} -o {output_tet};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_tig} -o {output_tig};
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_tob} -o {output_tob}
        """


rule combine_files:
    input:
        predictions_ami="CNN_predictions_ami.txt"
        predictions_amo="CNN_predictions_amo.txt"
        predictions_amp="CNN_predictions_amp.txt"
        predictions_az="CNN_predictions_az.txt"
        predictions_cefe="CNN_predictions_cefe.txt"
        predictions_ceft="CNN_predictions_ceft.txt"
        predictions_chlo="CNN_predictions_chlo.txt"
        predictions_cip="CNN_predictions_cip.txt"
        predictions_clin="CNN_predictions_clin.txt"
        predictions_col="CNN_predictions_col.txt"
        predictions_dor="CNN_predictions_dor.txt"
        predictions_ert="CNN_predictions_ert.txt"
        predictions_eryt="CNN_predictions_eryt.txt"
        predictions_fos="CNN_predictions_fos.txt"
        predictions_gen="CNN_predictions_gen.txt"
        predictions_imi="CNN_predictions_imi.txt"
        predictions_lev="CNN_predictions_lev.txt"
        predictions_mer="CNN_predictions_mer.txt"
        predictions_mox="CNN_predictions_mox.txt"
        predictions_nit="CNN_predictions_nit.txt"
        predictions_tet="CNN_predictions_tet.txt"
        predictions_tig="CNN_predictions_tig.txt"
        predictions_tob="CNN_predictions_tob.txt"
    output:
        final=temp("static/files/CNN/{genome}_FINAL.txt.CNN_predicitions.txt")
    run:
        shell("""echo "Antibiotic predicted_pheno" > {output.final}; awk '{{print "Amikacin", $0}}' {input.predictions_ami} >> {output.final}; awk '{{print "Amoxicillin", $0}}' {input.predictions_amo} >> {output.final}; awk '{{print "Ampicillin", $0}}' {input.predictions_amp} >> {output.final}; awk '{{print "Aztreonam", $0}}' {input.predictions_az} >> {output.final}; awk '{{print "Cefepime", $0}}' {input.predictions_cefe} >> {output.final}; awk '{{print "Ceftriaxone", $0}}' {input.predictions_ceft} >> {output.final}; awk '{{print "Chloramphenicol", $0}}' {input.predictions_chlo} >> {output.final}; awk '{{print "Ciprofloxacin", $0}}' {input.predictions_cip} >> {output.final}; awk '{{print "Clindamycin", $0}}' {input.predictions_clin} >> {output.final}; awk '{{print "Colistin", $0}}' {input.predictions_col} >> {output.final}; awk '{{print "Doripenem", $0}}' {input.predictions_dor} >> {output.final}; awk '{{print "Ertapenem", $0}}' {input.predictions_ert} >> {output.final}; awk '{{print "Erythromcyin", $0}}' {input.predictions_eryt} >> {output.final}; awk '{{print "Fosfomycin", $0}}' {input.predictions_fos} >> {output.final}; awk '{{print "Gentamicin", $0}}' {input.predictions_gen} >> {output.final}; awk '{{print "Imipenem, $0}}' {input.predictions_imi} >> {output.final}; awk '{{print "Levofloxacin", $0}}' {input.predictions_lev} >> {output.final}; awk '{{print "Meropenem", $0}}' {input.predictions_mer} >> {output.final}; awk '{{print "Moxifloxacin", $0}}' {input.predictions_mox} >> {output.final}; awk '{{print "Nitrofurantoin", $0}}' {input.predictions_nit} >> {output.final}; awk '{{print "Tetracycline", $0}}' {input.predictions_tet} >> {output.final}; awk '{{print "Tigecycline", $0}}' {input.predictions_tig} >> {output.final}; awk '{{print "Tobramycin", $0}}' {input.predictions_tob} >> {output.final}""")



rule DT_CNN:
    input:
        DT="static/files/decision_tree/{genome}_FINAL.txt.DT.predicitions.txt"
        CNN="static/files/CNN/{genome}_FINAL.txt.CNN_predicitions.txt"
    output:
        file="static/files/FINAL_predictions.txt"
    shell:
        """
        cat {input.DT} {input.CNN} > {output.file}
        """





rule send_results:
    input:
        file="static/files/FINAL_predictions.txt"
    params:
        name="static/files/decision_tree/{genome}.txt"
    output:
        result="static/files/results/{genome}.txt"
    shell:
        """
        mv {input.file} {params.name};
        echo "both models smk" >> {params.name};
        mv {params.name} {output.result}
        """
