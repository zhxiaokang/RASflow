# The main script to manage the subworkflows of RASflow

import yaml
import os

with open('configs/config_main.yaml') as yamlfile:
    config = yaml.load(yamlfile)

# Parameters to control the workflow

## Do you need to do quality control?
qc = config["QC"]
print("Is quality control required?\n", qc)

## Do you need to do trimming?
trim = config["TRIMMED"]
print("Is trimming requred?\n", trim)

## Which mapping reference do you want to use? Genome or transcriptome?
reference = config["REFERENCE"]
print("Which mapping reference will be used?\n", reference)

## Do you want to do Differential Expression Analysis (DEA)?
dea = config["DEA"]
print("Is DEA requred?\n", dea)

## Do you want to visualize the results of DEA?
visualize = config["VISUALIZE"]
print("Is visualization requred?\n", visualize)

# Double check with the user about the requested sub-workflows to be run
print("Please double check the information above\nDo you want to continue? (y/n)")
check_flow = input()
if check_flow == "y":
    pass
else:
    os._exit(0)

# Start the workflow
print("Start RASflow!")

if qc:
    # Double check that the user really wants to do QC instead of forgetting to change the param after doing QC
    print("Are you sure that you want to do Quality Control?\n If yes, type 'y'; if not, type 'n' and set 'QC' to 'no' in the config file")
    qc_2nd = input()
    if qc_2nd == "y":
        os.system("nice -5 snakemake -s workflow/quality_control.rules 2>&1 | tee logs/log_quality_control.txt")
        print("Quality control is done!\n Please check the report and decide whether trimming is needed\n Please remember to turn off the QC in the config file!")
        os._exit(0)
    else:
        os._exit(0)
else:
    if trim:
        os.system("nice -5 snakemake -s workflow/trim.rules 2>&1 | tee logs/log_trim.txt")
        print("Trimming is done!")
    else:
        print("Trimming is not required")

    print("Start mapping using ", reference, " as reference!")

    if reference == "transcriptome":
        os.system("nice -5 snakemake -s workflow/quantify_trans.rules 2>&1 | tee logs/log_quantify_trans.txt")
    elif reference == "genome":
        os.system("nice -5 snakemake -s workflow/align_count_genome.rules 2>&1 | tee logs/log_align_count_genome.txt")

    if dea:
        print("Start doing DEA!")
        if reference == "transcriptome":
            os.system("nice -5 snakemake -s workflow/dea_trans.rules 2>&1 | tee logs/log_dea_trans.txt")
        elif reference == "genome":
            os.system("nice -5 snakemake -s workflow/dea_genome.rules 2>&1 | tee logs/log_dea_genome.txt")
        print("DEA is done!")

        if visualize:
            # Visualization can only be done on gene-level
            if reference == "genome":
                pass
            elif reference == "transcriptome":
                gene_level = config["GENE_LEVEL"]
                if gene_level:
                    pass
                else:
                    print("Sorry! RASflow currently can only visualize on gene-level")
                    os._exit(1)

            print("Start visualization of DEA results!")
            os.system("nice -5 snakemake -s workflow/visualize.rules 2>&1 | tee logs/log_visualize.txt")
            print("Visualization is done!")
            print("RASflow is done!")
        else:
            print("Visualization is not required and RASflow is done!")
    else:
        print("DEA is not required and RASflow is done!")
