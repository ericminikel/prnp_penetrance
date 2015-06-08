#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

'''
Three inputs:
1. Genotypes VCF
2. Annotated sites VCF
3. Metadata file (.TSV)

Three outputs in specified directory:
1. samples.txt file for making IGV screenshots
2. sample_locus.txt file for making IGV screenshots
3. screenshot_annotations.txt file for annotating IGV screenshots
'''

import sys
import re
import vcf
import gzip
import argparse
import os.path
import re

parser = argparse.ArgumentParser(description='Process individuals with variants from a VCF.')
parser.add_argument('--genotypes', '-g', dest='genotypes', action='store',
                    help='Genotypes VCF file')
parser.add_argument('--annotations', '-a', dest='annotations', action='store',
                    help='Annotated sites VCF file')
parser.add_argument('--metadata', '-m', dest='metadata', action='store',
                    help='Metadata .tsv file')
parser.add_argument('--pass', dest='passonly', action='store_const',
                   const=True, default=False, help='Only PASS variants. Not recommended due to multi-allelic sites' )
parser.add_argument('--maxmaf', dest='maxmaf', action='store', type=float,
                   default=1.0, help='Maximum MAF to bother with' )
parser.add_argument('--mingq', dest='mingq', action='store', type=int,
                   default=0.0, help='Minimum GQ to bother with' )
parser.add_argument('--minad', dest='minad', action='store', type=int,
                   default=0.0, help='Minimum AD of alt allele' )
parser.add_argument('--minab', dest='minab', action='store', type=float,
                   default=0.0, help='Minimum allelic balance' )
parser.add_argument('--outdir', dest='outdir', action='store',
                   default='./', help='Output directory [default ./]' )
parser.add_argument('--region', dest='region', action='store',
                   default='', help='Genomic region of interest in CHROM:START-STOP format [default all]' )
parser.add_argument('--transcript', dest='transcript', action='store',
                   default='', help='Transcript of interest [default not tested]' )
args = parser.parse_args()

real_prnp_tx = "ENST00000379440" # of the many Ensembl transcripts for PRNP, this is the "real" one, that we care about
symbols = re.compile('[^\w]') # a pattern matcher for any symbols I don't want in filenames

def smartopen(path):
    '''
    Open a path with gzip if gzipped, regular open otherwise.
    '''
    if path[-3:] == ".gz":
        return gzip.open(path)
    else:
        return open(path)

def aa_change(amino_acids,protein_pos):
    if len(amino_acids) == 0: # non-coding variant
        return ''
    elif len(amino_acids) == 1: # synonymous variant e.g. "A" means alanine-to-alanine
        return amino_acids+protein_pos+amino_acids # e.g. "A117A"
    else:
        amino_acids = amino_acids.replace("*","X") # use X instead of * to indicate stop codons
        aas = amino_acids.split('/')
        return aas[0] + protein_pos + aas[1]

# open vcf readers for the two VCFs, and read all the metadata
annos = vcf.Reader(smartopen(args.annotations),filename='ignore',compressed=False,strict_whitespace=True)
genos = vcf.Reader(smartopen(args.genotypes),filename='ignore',compressed=False,strict_whitespace=True)
with open(args.metadata) as m:
    metaheader = m.readline()
    metalines = m.readlines()

# parse metadata into a dictionary of dictionaries.
colnames = metaheader.strip().split('\t')
meta = {}
for metaline in metalines:
    sub_dict = {}
    cols = metaline.strip().split('\t')
    for i in range(len(cols)):
        sub_dict[colnames[i]] = cols[i]
    meta[cols[0]] = sub_dict

# store output as lists for now, only write to files at end
sample_lines = []
sample_locus_lines = []
screenshot_annotation_lines = []

# add a header to the screenshot annotations
screenshot_annotation_lines.append("\t".join(["png_filename","allele_id","hgvsc","aa_code","sample","exact_gt_correct","comments",
                    "this_alt_ad","trans_allele_ad","this_alt_ab","this_gt_gq"]))

# strategy: iterate over annotations, and pull out the corresponding row from genotypes

# first, pre-load the annotation records since they are small
anno_records = list(annos)
csq_collist = annos.infos['CSQ'].desc.split()[-1] # extract the pipe-delimited list of CSQ column names from the annotation VCF
csq_colnames = csq_collist.split('|')

for anno_record in anno_records:
    geno_record = genos.next() # pull corresponding genotype record
    # check that the two are referring to the same variant(s)
    assert str(anno_record) == str(geno_record), "Annotation and genotype lines do not match:\n%s\n%s\n"%(str(anno_record),str(geno_record))
    # iterate over multiple alleles
    for alt in geno_record.ALT:
        this_alt_allele_index = geno_record.ALT.index(alt) # index of this particular allele in comma-separated INFO fields
        this_alt_allele_number = geno_record.ALT.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
        nominal_af = float(geno_record.INFO['AF'][this_alt_allele_index]) # extract allele freq for this allele
        # decide whether to process this allele
        if nominal_af > args.maxmaf:  # only continue for alleles in the desired AF range
            continue
        if args.passonly and len(geno_record.FILTER) > 0: # if this is a filtered site and user wants PASS-only
            continue # then move on (but this feature is not recommended, see above)
        for sample in geno_record.samples: # iterate over all individuals int he VCF
            # decide whether to process this sample
            if sample['GT'] is None: # no-calls apparently come through as None instead of ./.
                # if you call sample.gt_alleles on them, PyVCF tries to do None.split() and
                # throws an Attribute Error. so just ignore these.
                continue
            if this_alt_allele_number not in map(int,sample.gt_alleles): # if this sample does not have this allele
                continue # move on
            # if you get here, the sample DOES have the allele. now process accordingly.
            this_alt_ad = int(sample['AD'][map(int,sample.gt_alleles).index(this_alt_allele_number)]) # extract allelic depth for this allele
            if this_alt_ad < args.minad: # if it's below the minimum, ignore and move on
                print "SKIPPING", str(this_alt_ad), str(args.minad)
                continue
            this_gt_gq = int(sample['GQ']) # extract GQ for this sample's genotype
            if this_gt_gq < args.mingq: # if it's below the minimum, ignore and move on
                print "SKIPPING", str(this_gt_gq), str(args.mingq)
                continue
            this_gt_dp = int(sample['DP']) # extract total depth for this genotype
            this_alt_ab = float(this_alt_ad) / float(this_gt_dp) # calculate allelic balance for this allele
            if this_alt_ab < args.minab: # if it's below the minimum, ignore and move on
                print "SKIPPING", str(this_alt_ab), str(args.minab)
                continue
            trans_allele_ad = int(sample['AD'][1-map(int,sample.gt_alleles).index(this_alt_allele_number)]) # extract allelic depth for opposite allele
            # note this may not just be this_gt_dp - this_alt_ad because there may be a minority of reads supporting a third allele
            print sample.sample, this_alt_ad, this_gt_gq, this_gt_dp, this_alt_ab
            # if you make it to this point, we are interested in this allele.
            # so now pull out annotations.
            if not meta.has_key(sample.sample): # first check there is metadata for the sample
                sys.stderr.write("No metadata found for %s"%(sample.sample)) # if not, tell me which sample and print the details...
                sys.stderr.write("\t".join([str(geno_record), sample.sample, this_alt_ad, this_gt_gq, this_gt_dp, this_alt_ab]))
                continue # ... but then move on gracefully
            # contstruct allele id and file-friendly identifiers
            allele_id = geno_record.CHROM + ":" + str(geno_record.POS) + "_" + geno_record.REF + ">" + str(alt)
            allele_id_filename_friendly = geno_record.CHROM + "_" + str(geno_record.POS) + "_" + geno_record.REF + "_" + str(alt)
            sample_id_file_friendly = re.sub(symbols,'_',sample.sample)
            locus_range = geno_record.CHROM + ":" + str(geno_record.POS) + "+50"
            sample_locus_id = geno_record.CHROM + "_" + str(geno_record.POS) + "_" + sample_id_file_friendly + "_" + str(alt) + "_" + str(this_alt_ad)
            png_filename = sample_locus_id + ".png"
            exact_gt_correct = '' # empty field for user to fill in on screenshot annnotations
            comments = '' # empty field for user to fill in on screenshot annnotations
            # extract metadata
            bampath     = meta[sample.sample]['bam']
            projectid   = meta[sample.sample]['ProjectID']
            projectname = meta[sample.sample]['ProjectName']
            # extract annotations
            all_vep_annos = anno_record.INFO['CSQ'] # returns a list of annotations, one for each Ensembl transcript
            csq = None # initialize what will hopefully be a dictionary of VEP annotations
            for vep_anno in all_vep_annos: # loop to grab just the annotations for the transcript and alelle we care about
                csq = dict(zip(csq_colnames,vep_anno.split('|'))) # create a dict of the current entry
                # of all the VEP annotations, use only the ones corresponding to *this allele* on *our transcript of interest*
                if csq['Feature'] == real_prnp_tx and csq['ALLELE_NUM'] == str(this_alt_allele_number):
                    break # otherwise, keep looping until you run out of VEP annotations
            assert csq is not None, "No VEP annotations matched transcript %s and allele number %s\n%s"%(real_prnp_tx,str(this_alt_allele_number),str(anno_record))
            aa_code = aa_change(csq['Amino_acids'],csq['Protein_position']) # create a one-letter amino acid change code, e.g. "E200K"
            # extract other CSQ annotations of interest
            hgvsc = csq['HGVSc']
            # check if the bam exists.
            if not os.path.isfile(bampath):
                newbampath = re.sub('v[0-9]','current',bampath) # if it doesn't, see if there's a current symlink
                if os.path.isfile(newbampath): # if yes, just use that...
                    bampath = newbampath
                else: # if not, we have no BAM and can't process this sample
                    sys.stderr.write("BAM does not exist: %s\n"%(bampath))
                    continue # should discuss how best to record the fail rate here.
            # add data to lists
            sample_lines.append("\t".join([sample_id_file_friendly,bampath]))
            sample_locus_lines.append("\t".join([sample_locus_id,locus_range,sample_id_file_friendly]))
            screenshot_annotation_lines.append("\t".join(map(str,[png_filename,allele_id,hgvsc,aa_code,sample,exact_gt_correct,comments,
                projectid,projectname,this_alt_ad,trans_allele_ad,this_alt_ab,this_gt_gq])))
            print allele_id, aa_code, str(geno_record), sample.sample, this_alt_ad, this_gt_gq, this_gt_dp, this_alt_ab, bampath, projectid, projectname


with open(args.outdir+'samples.txt',mode='wb') as s:
    for line in sample_lines:
        s.write(line+"\n")

with open(args.outdir+'sample_locus.txt',mode='wb') as sl:
    for line in sample_locus_lines:
        sl.write(line+"\n")

with open(args.outdir+'screenshot_annotations.txt',mode='wb') as sa:
    for line in screenshot_annotation_lines:
        sa.write(line+"\n")


