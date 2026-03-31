from os import listdir
from os.path import isfile, join
from re import findall

configfile: 'config.yaml'

output_path = config['output_path'] #Path to the output folder
base_path = config['gwas_summstat_path'] #Path to the folder containing the formatted GWAS summary statistics for each trait

phenotype = [findall(r'(\w+)_gwas_ssf.h.tsv.gz', f)[0] for f in listdir(base_path) if isfile(join(base_path, f)) and findall(r'(\w+)_gwas_ssf.h.tsv.gz', f)] 
#Extracting the trait names from the formatted GWAS summary statistics files

#Writing the all rule which defines the final outputs of 1) hm3.sumstats.gzfile for each phenotype,2) SNP-heritability estimation, 3) Genetic correlation estimation with HCM disease state
rule all:
    input: 
        expand(output_path+'{p}_ldsc.hm3.sumstats.gz', p=phenotype),
        expand(output_path+'{p}_ldsc.hsq.log', p=phenotype), 
        expand(output_path+'{p}_ldsc.hcm.gcorr.log', p=phenotype)
        #,
        # output_path+'hsq_summary.tsv',
        # output_path+'gcorr_summary.tsv'

#Writing the rule to munge the GWAS summary statistics
rule munge:
    input: base_path+'{phenotype}_gwas_ssf.h.tsv.gz'
    output: output_path+'{phenotype}_ldsc.hm3.sumstats.gz'
    resources:
        mem_mb = 16000
    params:
        output_name= output_path+'{phenotype}_ldsc.hm3',
        merge_alleles_snplist = config['merge_alleles_snplist']
    shell:'''
    module load Miniforge3/24.1.2-0
    eval "$(conda shell.bash hook)"
    conda activate ldsc39

    /well/PROCARDIS/jchan/bin/ldsc/ldsc39/munge_sumstats.py \
    --sumstats {input} \
    --out {params.output_name} \
    --merge-alleles {params.merge_alleles_snplist}
    '''

#Writing the rule to estimate SNP-heritability
rule snp_heritability:
    input:output_path+'{phenotype}_ldsc.hm3.sumstats.gz'
    output: output_path+'{phenotype}_ldsc.hsq.log'
    resources:
        mem_mb = 16000
    params:
        output_name= output_path+'{phenotype}_ldsc.hsq',
        ld_path = config['ld_folder']
    shell:'''
        module load Miniforge3/24.1.2-0
        eval "$(conda shell.bash hook)"
        conda activate ldsc39

        /well/PROCARDIS/jchan/bin/ldsc/ldsc39/ldsc.py \
        --h2 {input} \
        --ref-ld-chr {params.ld_path} \
        --w-ld-chr {params.ld_path} \
        --out {params.output_name}
    '''

#Writing the rule to estimate genetic correlation with HCM disease state
rule genetic_correlation:
    input: output_path+'{phenotype}_ldsc.hm3.sumstats.gz'
    output: output_path+'{phenotype}_ldsc.hcm.gcorr.log'
    resources:
        mem_mb = 16000
    params:
        input_name= output_path+'{phenotype}_ldsc.hm3.sumstats.gz,'+config['tadros25_gwas_summstats'],
        output_name= output_path+'{phenotype}_ldsc.hcm.gcorr',
        ld_path=config['ld_folder']
    shell:'''
        module load Miniforge3/24.1.2-0
        eval "$(conda shell.bash hook)"
        conda activate ldsc39

        /well/PROCARDIS/jchan/bin/ldsc/ldsc39/ldsc.py \
        --rg {params.input_name} \
        --ref-ld-chr {params.ld_path} \
        --w-ld-chr /{params.ld_path} \
        --out {params.output_name}
        '''

# #Summarise the entire outputs
# rule ldsc_summariser:
#     input:
#         expand(output_path+'{p}_ldsc.hsq.log', p=phenotype), 
#         expand(output_path+'{p}_ldsc.hcm.gcorr.log', p=phenotype)
#     output:
#         hsq_summary=output_path+'hsq_summary.tsv',
#         gcorr_summary=output_path+'gcorr_summary.tsv'
#     params:
#         input_path= output_path
#     resources:
#         mem_mb=4000
#     shell:'''
#         module purge
#         module load R/4.2.2-foss-2022b
#         Rscript 1_LDSC_summariser.R {params.input_path}
#     '''
