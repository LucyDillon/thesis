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
         cnn_model="Training_models/gentamicin_egg_model.h5"
    output:
          "CNN_predictions.txt"
    shell:
        """
        python3 eggnog2cnnPredictions.py -file {input.eggnog_output} -model {input.cnn_model} -o {output}
        """

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
