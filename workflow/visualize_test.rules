import os

configfile: "configs/config_visualize.yaml"

dea_path = config["DEAPATH"]
norm_path = config["NORMPATH"]
control = config["CONTROL"][0]
treat = config["TREAT"][0]
output_path = config["OUTPATH"]

rule end:
    input:
        volcano = output_path + "/volcano_plot_" + control + "_" + treat + ".pdf",
        heatmap = output_path + "/heatmap_" + control + "_" + treat + ".pdf"

rule plot:
    input:
        dea = dea_path + "/dea_" + control + "_" + treat + ".tsv",
        norm_control = norm_path + "/" + control + "_norm.tsv",
        norm_treat = norm_path + "/" + treat + "_norm.tsv"
    output:
        volcano = output_path + "/volcano_plot_" + control + "_" + treat + ".pdf",
        heatmap = output_path + "/heatmap_" + control + "_" + treat + ".pdf"
    shell:
        "Rscript scripts/visualize.R"

