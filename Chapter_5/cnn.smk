singularity: "./amrmlpipeline_6.sif"
GENOMES, = glob_wildcards("static/files/{genome}.fna")

rule all:
    input:
        expand("static/files/results/{genome}.txt", genome=GENOMES)

rule prodigal:
    input:
        "static/files/{genome}.fna"
    output:
        temp("static/files/prodigal_output/{genome}.faa")
    shell:
        '''
        prodigal -i {input} -a {output} -f gff
        '''

rule diamond:
    input:
        faa=temp("static/files/prodigal_output/{genome}.faa"),  # Mark as temporary
        db="data/eggnog_proteins.dmnd"
    output:
        temp("static/files/diamond_annotation/{genome}_Diamond.faa")
    shell:
        '''
        diamond blastp -d {input.db} -q {input.faa} --more-sensitive --threads 2 -e 0.001000 -o {output} --top 3 --query-cover 0 --subject-cover 0
        '''

rule Edit_diamond:
    input:
        "static/files/diamond_annotation/{genome}_Diamond.faa"
    output:
        "static/files/diamond_annotation/{genome}_DiamondReduced.faa"
    run:
        shell("""cat {input} | awk 'BEGIN{{prev = ""}};{{if ($1 != prev) {{ prev=$1; print $1"\t"$2"\t"$11"\t"$12 }} }}' > {output}""")

rule eggnog:
    input:
	faa="static/files/diamond_annotation/{genome}_DiamondReduced.faa",
        db="data/"
    output:
	"{genome}.emapper.annotations"  # Correct the output file name here
    params:
	prefix="{genome}"
    shell:
	"""
	emapper.py --data_dir {input.db} --annotate_hits_table {input.faa}  --no_file_comments --override --output {params.prefix} --cpu 8
        """


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
        final=temp("static/files/decision_tree/{genome}_FINAL.txt.CNN_predicitions.txt")
    run:
        shell("""echo "Antibiotic predicted_pheno" > {output.final}; awk '{{print "Amikacin", $0}}' {input.predictions_ami} >> {output.final}; awk '{{print "Amoxicillin", $0}}' {input.predictions_amo} >> {output.final}; awk '{{print "Ampicillin", $0}}' {input.predictions_amp} >> {output.final}; awk '{{print "Aztreonam", $0}}' {input.predictions_az} >> {output.final}; awk '{{print "Cefepime", $0}}' {input.predictions_cefe} >> {output.final}; awk '{{print "Ceftriaxone", $0}}' {input.predictions_ceft} >> {output.final}; awk '{{print "Chloramphenicol", $0}}' {input.predictions_chlo} >> {output.final}; awk '{{print "Ciprofloxacin", $0}}' {input.predictions_cip} >> {output.final}; awk '{{print "Clindamycin", $0}}' {input.predictions_clin} >> {output.final}; awk '{{print "Colistin", $0}}' {input.predictions_col} >> {output.final}; awk '{{print "Doripenem", $0}}' {input.predictions_dor} >> {output.final}; awk '{{print "Ertapenem", $0}}' {input.predictions_ert} >> {output.final}; awk '{{print "Erythromcyin", $0}}' {input.predictions_eryt} >> {output.final}; awk '{{print "Fosfomycin", $0}}' {input.predictions_fos} >> {output.final}; awk '{{print "Gentamicin", $0}}' {input.predictions_gen} >> {output.final}; awk '{{print "Imipenem, $0}}' {input.predictions_imi} >> {output.final}; awk '{{print "Levofloxacin", $0}}' {input.predictions_lev} >> {output.final}; awk '{{print "Meropenem", $0}}' {input.predictions_mer} >> {output.final}; awk '{{print "Moxifloxacin", $0}}' {input.predictions_mox} >> {output.final}; awk '{{print "Nitrofurantoin", $0}}' {input.predictions_nit} >> {output.final}; awk '{{print "Tetracycline", $0}}' {input.predictions_tet} >> {output.final}; awk '{{print "Tigecycline", $0}}' {input.predictions_tig} >> {output.final}; awk '{{print "Tobramycin", $0}}' {input.predictions_tob} >> {output.final}""")

rule send_results:
    input:
        file="static/CNN_predictions.txt"  
    params:
        name="static/"
    output:
        result="static/files/results/{genome}_CNN.txt"
    shell:
        """
        mv {input.file} {params.name};
        echo "CNN model smk" >> {params.name};
        mv {params.name} {output.result}
        """
