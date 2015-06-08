#!/usr/bin/env python

import os
import sys
import re
from PIL import Image

# paths
annotations = 'data_nosync/screenshots_annotated.tsv'
screenshot_dir = 'data_nosync/igv'
destination_dir = 'supplement/igv'
metadata = 'data_nosync/ExAC_meta_Final.tsv'

# constants
crop_left = 200
crop_top = 0
crop_right = 40
crop_bottom = 110

# parse the ExAC metadata file
meta = {}
with open(metadata, mode='r') as m:
    colnames = m.readline().strip('\n').lower().split('\t')
    for line in m.readlines():
        data = dict(zip(colnames, line.strip('\n').split('\t')))
        meta[data['vcf_sampleid']] = data

# crop and save images with new filenames
succeeded = 0
skipped = 0
failed = 0
with open(annotations, mode='r') as a:
    colnames = a.readline().strip('\n').lower().split('\t')
    dup_counter = 1
    last_allele_id = None
    for line in a.readlines():
        data = dict(zip(colnames, line.strip('\n').split('\t')))
        sampleid = re.sub(",.*","",re.sub("Call\(sample=","",data['details']))
        include = meta.get(sampleid, {'include': 'NO'}).get('include') == "YES"
        if data['exact_gt_correct'] != 'y' or not include: # only create cropped screenshots of the variant calls deemed correct and included in v0.3 release
            skipped += 1
            continue
        image_filename = os.path.join(screenshot_dir, data['png_filename'])
        allele_id = re.sub(r'[^\w]', '_', data['allele_id']) # replace CHROM:POS_REF>ALT with CHROM_POS_REF_ALT
        aa_code = data['aa_code']
        if allele_id == last_allele_id:
            dup_counter += 1
        else:
            dup_counter = 1
        if os.path.isfile(image_filename):
            img = Image.open(image_filename)
            width, height = img.size
            cropped = img.crop((crop_left, crop_top, width-crop_right, height-crop_bottom))
            destination_filename = os.path.join(destination_dir, allele_id + "_" + aa_code + "_" + str(dup_counter) + ".png")
            cropped.save(destination_filename)
            succeeded += 1
        else:
            failed += 1
            sys.stderr.write('Failed on this entry: ' + str(data) + '\n')
        last_allele_id = allele_id

print "Successfully cropped %s images, skipped %s images, failed on %s images."%(succeeded, skipped, failed)


