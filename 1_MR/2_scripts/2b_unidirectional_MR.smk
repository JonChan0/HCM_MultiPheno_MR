''' 
Snakemake script to run Bidirectional Mendelian Randomisation in parallel over multiple traits.
Author: Jonathan Chan
Date: 2025-02-20
'''

configfile: 'config/config_unidirectionalMR_hcmrMR.yaml'

#Primary trait and secondary trait are like trait 1 and trait 2 - they take turns being exposure and outcome respectively
rule all:
    input:
        # Direction 2: Primary Trait -> Secondary Trait
        expand(config['output_path'] + '{exposure}_to_{outcome}/{exposure}_mr_results_mainline.tsv',
               exposure=config['primary_trait'], outcome=config['secondary_traits']),
        expand(config['output_path'] + '{exposure}_to_{outcome}/{exposure}_mr_results_sensitivity_oddsratio.tsv',
               exposure=config['primary_trait'], outcome=config['secondary_traits'])


# A single, generalised rule to run MR in either direction.
# Wildcards {exposure} and {outcome} are used to fetch the correct data from the config.
rule run_mr:
    input:
        # Use a lambda function to dynamically find the correct path from the config file.
        exposure_path = lambda wildcards: config['trait_info'][wildcards.exposure]['path'],
        outcome_path = lambda wildcards: config['trait_info'][wildcards.outcome]['path']
    output:
        # Output paths are generated using the wildcards.
        tsv_output = config['output_path'] + '{exposure}_to_{outcome}/{exposure}_mr_results_mainline.tsv',
        tsv2_output = config['output_path'] + '{exposure}_to_{outcome}/{exposure}_mr_results_sensitivity_oddsratio.tsv'
    params:
        output_path = config['output_path'],
        exposure = '{exposure}',
        outcome = '{outcome}',
        # This lambda function handles the complex logic for selecting the right instruments.
        # If the exposure is the primary trait, it looks up the instruments based on the outcome.
        # Otherwise, it uses the standard instrument list for that exposure.
        selected_instruments = lambda wildcards: config['trait_info'][wildcards.exposure]['instruments_as_exposure'],
        ld_clumping = lambda wildcards: config['trait_info'][wildcards.exposure].get('ld_clump', 'TRUE'),
        ld_eur_bed_file = config['ld_eur_bed_file'],
        plink_binary_path = config['plink_binary_path']
    resources:
        mem_mb = 32000
    shell:
        '''
        module purge
        module use -a /apps/eb/2022b/skylake/modules/all
        module load R/4.2.2-foss-2022b
        Rscript 2_TwoSampleMR.R {input.exposure_path} {input.outcome_path} {params.selected_instruments} {params.output_path} {params.exposure} {params.outcome} {params.ld_clumping} {params.ld_eur_bed_file} {params.plink_binary_path}
        '''
