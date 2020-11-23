from wgs.utils import helpers as h
import pandas as pd
import pyranges as pr
import statistics
from scipy.stats import norm
import sys

def label(cn, mu):
'''
label a copy number state value based on fitted normal distribution mu
'''
     if cn < 1:
         return 'HDEL'
     elif 1 <= cn <= mu-1:
         return 'DEL'
     elif mu+1 <= cn < 6:
         return 'AMP'
    elif cn >= 6:
        return 'HLAMP'
    else:
        return "NO_CHANGE"


def label_cn_states(cn):
    '''
    label array of copy number states cn
    '''
    mu, _ = norm.fit(cn.state.tolist())
    cn["state"] = cn.state.apply(lambda cn: label(cn, mu))
    return cn

def anno_maf_with_hmmcopy(maf, hmmcopy, cna_annotation, genes=None):
    '''
    take intersection of maf snvs and hmmcopy and label matches
    '''
    genes = pd.read_csv(maf, sep="\t", usecols=["Hugo_Symbol", "Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position"])
    genes = genes.rename(columns={"Start_Position": "Start", "End_Position":"End"})
    genes = genes[genes.Hugo_Symbol != "Unknown"]

    genes = pr.PyRanges(genes)

    hmmcopy = pd.read_csv(hmmcopy, chunksize=100000, usecols=["chr", "start", "end", "state"])

    output = []

    for chunk in hmmcopy:
        chunk = chunk.rename(columns={"chr": "Chromosome", "start":"Start", "end": "End"})
        chunk = chunk.groupby(["Chromosome", "Start", "End"])
        chunk = chunk['state']
        chunk = chunk.aggregate({'state': statistics.median})
        chunk = chunk.reset_index()
        # chunk = chunk[chunk.state != 2.0]

        chunk = pr.PyRanges(chunk)
        overlaps = genes.join(chunk).as_df()
        output.append(overlaps)

    final = pd.concat(output)
    final = label_cn_states(final)
    # print(final)
    final = final.drop_duplicates(['Hugo_Symbol','Tumor_Sample_Barcode'],keep= 'last')
    final.to_csv(cna_annotation, sep="\t", index=False)


anno_maf_with_hmmcopy(sys.argv[1], sys.argv[2], sys.argv[3])
# anno_maf_with_hmmcopy("/juno/work/shah/abramsd/CODE/wgs_cohort_qc_TESTRUNS/MAFS/sample_SA1090/strelka_indel.maf", 
#     "/work/shah/tantalus/SC-2640/results/hmmcopy/A96213A_reads.csv.gz", "/juno/work/shah/abramsd/CODE/wgs_cohort_qc_TESTRUNS/anno.tsv")
