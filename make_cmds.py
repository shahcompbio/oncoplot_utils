
import os
import pandas as pd
import yaml
from wgs.utils import helpers

def vcf2maf(yaml, cmds_file):
    out_dir = "/juno/work/shah/abramsd/dlp_cohort_oncoplot/DATA/MAFS_COMPLETE"
    yaml = helpers.load_yaml_flat(yaml)
    museq = {label: data["museq_vcf"] for label, data in yaml}
    strelka = {label: data["strelka_vcf"] for label, data in yaml}
    strelka_indel = {label: data["strelka_indel_vcf"] for label, data in yaml}

    museq = pd.DataFrame.from_dict(museq, columns=["label", "vcf"]).transpose()
    strelka = pd.DataFrame.from_dict(strelka, columns=["label", "vcf"]).transpose()
    strelka_indel = pd.DataFrame.from_dict(strelka_indel, columns=["label", "vcf"]).transpose()
    
    data = pd.concat([museq, strelka, strelka_indel])
    data["cmd"] = ["mkdir {}; bsub -e e0 -W 96:00 -M 100G bash -c '  vcf2maf.pl --input-vcf {} --output-maf {}--ref-fasta /juno/work/shah/reference/genomes/GRCh37-lite/GRCh37-lite.fa --filter-vcf 0 --vep-path /home/abramsd/miniconda3/envs/r-environment/bin --buffer-size 10 --vep-data /work/shah/reference/vep"] * len(data.index)
    data["maf"] =  data.apply(lambda r: os.path.join(out_dir, r.label, r.vcf.replace(".vcf.gz", ".maf")), axis=1)
    data["dir"] =  data.apply(lambda r: os.path.join(out_dir, r.label), axis=1)

    data["cmd"] = data.apply(lambda r: r.cmd.format(r.dir, r.vcf, r.vcf), axis=1)    
    
    data.cmd.to_csv(cmds_file, index=False)
                    
        
def merge_library_mafs(mafs, cmds_file):
    out_dir = "/juno/work/shah/abramsd/dlp_cohort_oncoplot/DATA/MAFS_COMPLETE"
    mafs=pd.read_csv(mafs, sep="\t", names=["maf"])
    mafs["label"] = mafs.maf.apply(lambda m: m.split("/")[8])
    cmd = "bsub bash -c 'python  /juno/work/shah/abramsd/CODE/wgs/wgs/workflows/cohort_qc/combine.py  \"{}\" \"{}\" {}"
    for group in mafs.groupby("label"):
        merged = os.path.join(out_dir, label, "combined.maf")
        cmd = cmd.format(group.mafs.tolist(), group.labels.tolist() * 3, merged)
        with open(cmds_file) as f:
          f.write(cmd + "\n")
          
          
def annotate(library_mafs, hmmcopy_reads, cmds_file):
   
    a = pd.read_csv(library_mafs, names=["maf"])
    
    a["cmd"] = [" bsub bash -c 'python /juno/work/shah/abramsd/CODE/wgs/wgs/workflows/cohort_qc/annotate_hmm_with_genes.py {} {} {}'"] * len(a.index)

    a["label"] = a.maf.apply(lambda m: m.split("/")[8])

    a["anno"] =  a.apply(lambda r: os.path.join("/juno/work/shah/abramsd/RESULTS/dlp_cohort_oncoplot/cna_annotation", r.label, "cna_annotation.tsv"), axis=1)

    a["dir"] =  a.apply(lambda r: os.path.join("/juno/work/shah/abramsd/RESULTS/dlp_cohort_oncoplot/cna_annotation", r.label), axis=1)

    a["reads"] = a.label.apply(lambda l: reads[l])

    a["cmd"] = a.apply(lambda r: r.cmd.format( r.maf, r.reads, r.anno), axis=1)

    a.cmd.to_csv(cmds_file, sep="\t", index=False)

    
def merge_library_data(data, cmd_file, out_dir):
    data = pd.read_csv(data, names=["file"])
    data["label"] = data.file.apply(lambda f: f.split("/")[8])
    cmd = "bsub bash -c 'python  /juno/work/shah/abramsd/CODE/wgs/wgs/workflows/cohort_qc/combine.py  \"{}\" \"{}\" {}"
    files = data.file.tolist()
    labels = data.label.tolist()
    output = os.path.join(out_dir, "merged")
    
   
