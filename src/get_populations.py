#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.7.1-sqlite3-rtrees/bin/python

# example usage: 
# get_populations.py --pcs exac_all_pca.csv --weights pc_weights.txt --ped_1kg_path integrated_call_samples.20130502.ALL.ped --samples 60.5k_nospace.list > pops_60.5k.txt
import argparse

parser = argparse.ArgumentParser(description='Generate a null distribution of diversity scores')
parser.add_argument('--pcs', '-p', dest='pcs', action='store',
                    help='Principal components (path)', type=str)
parser.add_argument('--weights', '-w', dest='wts', action='store',
                    help='Weights (path)', type=str)
parser.add_argument('--ped_1kg_path', dest='ped_1kg_path', action='store',
                    help='Path to 1000 Genomes ped file', type=str)
parser.add_argument('--samples', dest='samples', action='store',
                    help='Path to list of sample IDs, one per line', type=str)
args = parser.parse_args()

def get_1kg_pops(ped_1kg_path):
    # ped_1kg_path: ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/integrated_call_samples.20130502.ALL.ped
    # population definitions: http://www.1000genomes.org/category/frequently-asked-questions/population
    popd = {} # dictionary where keys are IIDs and values are POPs
    with open(ped_1kg_path) as f:
        header = f.readline() # discard the PED file header
        for line in f.readlines():
            fields = line.strip().split('\t')
            popd[fields[1]] = fields[6] # IID = fields[1], POP = fields[6]
    return popd

def get_1kg_centroids(popd, pcs, n=9):
    pc_sums = {} # compute running sums and counts
    pc_counts = {} 
    for iid in pcs.keys(): # iterate over individual IDs
        if popd.has_key(iid):
            if not pc_sums.has_key(popd[iid]): # initialize upon first encountering a new population
                pc_sums[popd[iid]] = [0] * n # fill the sum of PCs with zeroes
                pc_counts[popd[iid]] = 0
            pc_sums[popd[iid]] = [x + y for x, y in zip(pc_sums[popd[iid]], pcs[iid])] # add this indiv's PCs to the running total
            pc_counts[popd[iid]] += 1
    pc_centroids = {} # now we'll take sums/counts to get means
    for pop in pc_sums.keys():
        pc_centroids[pop] = [pc_sum/pc_counts[pop] for pc_sum in pc_sums[pop]]
    return pc_centroids

def assign_1kg_centroids(pcs, pc_centroids, weights):
    nearest_pop = {} # dictionary where keys will be IIDs and values will be nearest population
    for iid in pcs.keys():
        min_dist = float("inf") # we'll iterate over the pops to see which centroid has the least distance to this person
        argmin_dist = '' # population for which distance is minimized
        for pop in pc_centroids.keys():
            curr_dist = euclid_dist(pcs[iid],pc_centroids[pop],weights)
            if curr_dist < min_dist:
                argmin_dist = pop
                min_dist = curr_dist
        nearest_pop[iid] = argmin_dist
    return nearest_pop

def summarize_pops(nearest_pop):
    pop_counts = {} # dictionary with POP as keys and count(distinct IID) as values
    for iid in nearest_pop.keys():
        if not pop_counts.has_key(nearest_pop[iid]):
            pop_counts[nearest_pop[iid]] = 0
        pop_counts[nearest_pop[iid]] += 1
    return pop_counts

def printdict(d,alpha=True):
    if alpha:
        for key in sorted(d.keys()):
            print str(key) + '\t' + str(d[key])
    else:
        for x, y in d.items():
            print str(x) + '\t' + str(y)

def get_pop_ac_and_an(vcfpath,vcf_header,popdict,chr,pos,ref,alt):
    '''
    Accepts a path to a (bgzipped, tabix-indexed) VCF file, a dictionary
    mapping individuals to populations, and chr,pos,ref,alt for one allele. 
    Looks at the VCF and returns a dict with the AC and AN for each population.
    '''
    vcf_line_string = get_vcf_line(vcfpath,chr,pos)
    pseudo_vcf_file = StringIO(vcf_header+vcf_line_string)
    vcf_reader = vcf.Reader(pseudo_vcf_file,'r')
    records = list(vcf_reader)
    assert len(records) > 0, "No records found for that allele."
    assert len(records) < 2, "VCF contains >1 record at that position."
    record = records[0] # now knowing there is exactly 1 record, take it.
    assert chr == record.CHROM, "Extracted contig name %s does not match input %s"%(record.CHROM,chr)
    assert pos == record.POS, "Extracted position %s does not match input %s"%(record.POS,pos)
    assert ref == record.REF, "Extracted REF allele %s does not match input %s"%(record.REF,ref)
    assert alt in record.ALT, "Specified ALT allele %s not found at this site. Alleles are: %s"%(alt,str(record.ALT))
    this_alt_allele_index = record.ALT.index(alt) # index of this particular allele in comma-separated INFO fields
    this_alt_allele_number = record.ALT.index(alt) + 1 # for GT fields, need allele number: 1, 2, etc. remember REF allele is 0.
    # initialize dictionaries to hold ac and an for each pop
    pop_ac = dict((key,0) for key in list(set(popdict.values())))
    pop_an = dict((key,0) for key in list(set(popdict.values())))
    for sample in record.samples:
        if sample['GT'] is None: # In PyVCF 0.6.4 no-calls come through as None instead of ./.
            # if you call sample.gt_alleles on them, PyVCF tries to do None.split() and
            # throws an Attribute Error. so just ignore these. (In PyVCF 0.6.7 this is fixed)
            continue
        else:
            iid = sample.sample.replace(" ","_") # convert space to underscore in IIDs as they appear in the PCA file
            if popdict.has_key(iid):
                pop_an[popdict[iid]] += 2 # AN for this indiv's population is +2 b/c two alleles called
                pop_ac[popdict[iid]]+= map(int,sample.gt_alleles).count(this_alt_allele_number)
            else:
                sys.stderr.write("FYI, %s is not in the population dictionary"%iid)
    return pop_ac, pop_an

def read_pcs(path,n=9):
    '''
    Read principal components from a CSV file at the specified path.
    First column is sample id, next n are principal components. Additional
    columns may be present but will be ignored.
    '''
    pcs = {}
    with open(path) as inf:
        header = inf.readline().strip().split(',')
        for line in inf.readlines():
            cols = line.strip().split(',')
            sampleid = cols[0]
            samplepcs = [float(col) for col in cols[1:(n+1)]]
            pcs[sampleid] = samplepcs
    return pcs

def read_weights(weightpath):
    '''
    Reads a list of weights (presumably PC eigenvalues) from
    a file, whitespace and/or newline separated.
    '''
    with open(weightpath) as f:
        filecontents = f.read() # gulp whole file
        weights = filecontents.split() # on any whitespace including \n, \t or ' '
        return map(float,weights) # convert all to numerics


def euclid_dist(coords1,coords2,weights=None):
    '''
    Given two equal-length lists of coordinates in multi-dimensional space,
    return the Euclidean distance between the two points.
    '''
    assert len(coords1) == len(coords2), "Coordinate vectors differ in length"
    squared_diffs = [(coords1[i] - coords2[i])**2 for i in range(0,len(coords1))]
    if weights is not None:
        assert len(weights) == len(squared_diffs), "Weight vector is different length than coordinate vectors"
        squared_diffs = [weights[i]*squared_diffs[i] for i in range(0,len(weights))]
    euclidean_distance = sum(squared_diffs)**.5
    return (euclidean_distance)

def mean_euclid_dist(samples,pcs,weights=None,warn=True):
    '''
    Given a list of samples and a dictionary of principal components, calculate
    the mean Euclidean distance between pairs of samples. For n samples this
    requires n choose 2 comparisons. If warn=False, then silently skip individuals
    absent from the pcs file.
    '''
    n = len(samples)
    assert n <= 2500, "Too computationally intensive to do more than 2500 samples"
    assert n > 1, "Mean distance not defined for <2 points"
    n_pairs = comb(n,2,exact=True)
    mean_euclid_dist = 0.0
    valid_pairs = 0
    # check if all samples are in the PCs file
    valid_samples = []
    for i in range(len(samples)):
        if pcs.has_key(samples[i]):
            valid_samples.append(samples[i])
        else:
            if warn:
                sys.stderr.write("Warning: sample ID \"%s\" not found in PCs file. Ignoring.\n"%(samples[i]))
    for pair in combinations(valid_samples,2):
        # average as you go
        mean_euclid_dist += euclid_dist(pcs[pair[0]],pcs[pair[1]],weights) / n_pairs
        valid_pairs += 1
    assert valid_pairs >= 1, "No valid pairs found"
    return mean_euclid_dist

if __name__ == '__main__':
    # read in the sample ids
    sampleids = []
    with open(args.samples) as f:
        for line in f.readlines():
            sampleids.append(line.strip("\n"))
    
    pcs = read_pcs(args.pcs)
    wts = read_weights(args.wts)
    popd = get_1kg_pops(args.ped_1kg_path)
    centroids = get_1kg_centroids(popd,pcs)
    nearest_pop = assign_1kg_centroids(pcs, centroids, wts)
    
    these_pops = dict((sampleid, nearest_pop[sampleid]) for sampleid in sampleids)
    
    for sampleid in these_pops.keys():
        print "\t".join([sampleid,these_pops[sampleid]])


