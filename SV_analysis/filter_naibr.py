###
# This code is to filter bedpe file generated from Naibr and make it
# compatible with SURVIVOR, with additional filtration of lengths and quality.
###

import argparse, os
import pandas as pd
#from natsort import natsorted

def add_arguments():
    parser = argparse.ArgumentParser(
    description='''This code is to filter bedpe file generated from Naibr
     and make it compatible with SURVIVOR, with additional filtration
     of lengths and quality.''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( "-f", "--file_name", default="test",
        help="bedpe file generated from NAIBR")
    parser.add_argument("-t", "--svtype", default="DEL",
        help="SV type to keep [DEL, DUP, INV]")
    parser.add_argument("-M", "--maxlength", default=10000000,type=int,
        help="Max SV length to keep")
    parser.add_argument("-m", "--minlength", default=10000,type=int,
        help="Min SV length to keep")
    parser.add_argument("-q", "--quantile", default=0.1, type=float,
        help="Remove SV with QC less than the quantile")
    parser.add_argument("-d", "--max_distance", default=10000,type=int,
        help="Max distance between SV to be merged")
    args = parser.parse_args()
    return (parser, args)

def naibr_to_linkedSV(df):
    types = df["info"].str.split(";").str[0].str.replace("Type=","")
    col12 = df["info"].str.replace("^.*?;","")
    lengths = abs(df.iloc[:,4]-df.iloc[:,1])
    df2 = pd.concat([df.iloc[:,0:6],types,df.iloc[:,6],lengths, df.iloc[:,7],df.iloc[:,10],col12], axis=1)
    df2.columns = pd.array(["chrom1", "start1", "stop1", "chrom2", "start2", "stop2","sv_type",
                            "sv_id", "sv_length", "qual_score", "filter", "info"])
    return(df2)

### The following 2 functions are only for deletions, inversion and duplication,
### with an assumption that start1 = end1, start2 = end2, chr1 = chr2
def bedpe_to_bed(df):
    df2 = pd.concat([df.iloc[:,0:2], df.iloc[:,5], df.iloc[:,6:]], axis=1)
    return(df2)

def bed_to_bedpe(df):
    df2 = pd.concat([df.iloc[:,0:2], df.iloc[:,1],df.iloc[:,0], df.iloc[:,2],df.iloc[:,2],
                     df.iloc[:,3:]], axis=1)
    col_names = df2.columns.tolist()
    col_names[0:6] = ["chrom1", "start1", "stop1", "chrom2", "start2", "stop2"]
    df2.columns = col_names
    return(df2)

def sort_natural(df):
    x = list(range(1, 23))
    x.extend(["X", "Y", "EBV"])
    x = ["chr" + str(c) for c in x]
    #natsorted(df.iloc[:,0].unique())) can be used instead of x
    df.iloc[:,0] = pd.Categorical(df.iloc[:,0], ordered=True, categories=x)
    df = df.sort_values(by=df.columns[0:2].tolist())
    return(df)


def main():
    parser, aa = add_arguments()
    file_name = aa.file_name.replace(".bedpe","")

    names = ["chrom1", "start1", "stop1", "chrom2", "start2", "stop2", "callID", "score", "a" , 'b', 'c','info']
    #For testing:
    #df = pd.read_csv('201026.BLR.M10_merge.ema_final.naibr_sv_calls.bedpe', sep='\t', header=0, names=names)
    df = pd.read_csv(aa.file_name, sep='\t', header=0, names=names)

    # START
    df = naibr_to_linkedSV(df)
    df = sort_natural(df)
    df = df[df.sv_type == aa.svtype]
    df = df[df.sv_length <= aa.maxlength]
    df = df[df.sv_length >= aa.minlength]
    df = df[df.qual_score >= df.qual_score.quantile(q=aa.quantile)]
    df = bedpe_to_bed(df)

    bedfile = file_name + "_filtered_sorted.bed"
    bedfile_merged = file_name + "_filtered_sorted_merged.bed"
    bedpefile_merged = file_name + "_filtered_sorted_merged.bedpe"
    vcffile_merged = file_name + "_filtered_sorted_merged.vcf"

    df.to_csv(bedfile, index=None, header=None, sep='\t')
    os.system("bedtools merge -i {} -c 4,5,6,7,8,9 -o first,collapse,max,mean,first,count -d {} > tmp_{}".format(bedfile,aa.max_distance,bedfile_merged))
    os.system("bedtools intersect -v -a tmp_{0} -b hg38_black_list.bed -b hg38_gap.bed > {0}".format(bedfile_merged))

    df = pd.read_csv(bedfile_merged,header=None, sep='\t')
    df = bed_to_bedpe(df)
    df.to_csv(bedpefile_merged, index=False, header=False, sep='\t')

    os.system("SURVIVOR bedpetovcf {} tmp_{}".format(bedpefile_merged, vcffile_merged))
    os.system("less tmp_{0} | sed s'/STRANDS=[0-9][- 0-9];//'g > {0}".format(vcffile_merged))
    os.system("rm {0} {1} tmp_{1} tmp_{2}".format(bedfile, bedfile_merged, vcffile_merged))

if __name__ == '__main__':
    main()
