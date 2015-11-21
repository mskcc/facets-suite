#!/opt/common/CentOS_6-dev/bin/current/python
descr = "Concatenate FACETS out files into tabular per-sample format and cncf files into a long table. Output is OUT and CNCF."
import argparse, subprocess, os, re, glob

parser = argparse.ArgumentParser(description = descr, formatter_class = argparse.RawTextHelpFormatter)
parser.add_argument('-c', '--cncfFiles', help = 'List cncf files', nargs = '+', required = False)
parser.add_argument('-o', '--outFiles', help = 'List out files', nargs = '+', required = False)
parser.add_argument('-p', '--outPrefix', help = 'Prefix output', required = False)
args = parser.parse_args()

### If infiles are not specified, check directory
if args.cncfFiles is None:
    cncfFiles = glob.glob('*cncf.txt') 
else:
    cncfFiles = args.cncfFiles

if args.outFiles is None:
    outFiles = glob.glob('*out')
else:
    outFiles = args.outFiles

### Check if outPrefix specified
if args.outPrefix is None:
    CNCF = '%s_CNCF.txt' % os.path.basename(os.getcwd())
    OUT = '%s_OUT.txt' % os.path.basename(os.getcwd())
else:
    CNCF = '%s_CNCF.txt' % args.outPrefix
    OUT = '%s_OUT.txt' % args.outPrefix

### Concatenate cncf files into CNCF.txt
concatCall = '(cat %s | head -1; cat *cncf.txt | egrep -v "^ID") > %s' % (' '.join(cncfFiles), CNCF)
subprocess.call(concatCall, shell = True, stdout = open(os.devnull, 'w'), stderr = subprocess.STDOUT)

### Read and parse out files
### Write to OUT.txt
header = 'Sample\tFacets\tsnp.nbhd\tndepth\tpurity_cval\tcval\tmin.nhet\tgenome\tPurity\tPloidy\tdipLogR\tloglik\tflags\n'
with open(OUT, 'w') as output:

    # Parse out files
    output.write(header)
    for out in outFiles:
        f = open(out, 'r')
        content = [x.strip('\n') for x in f.readlines()]
        contentDict = dict(item.split("=") for item in content if '=' in item)
        sample = [value.strip() for key, value in contentDict.items() if re.search('TAG', key)][0]
        facets = [value.strip() for key, value in contentDict.items() if re.search('Facets version', key)][0]
        snpnbhd = [value.strip() for key, value in contentDict.items() if re.search('snp.nbhd', key)][0]
        ndepth = [value.strip() for key, value in contentDict.items() if re.search('ndepth', key)][0]
        puritycval = [value.strip() for key, value in contentDict.items() if re.search('purity_cval', key)][0]
        cval = [value.strip() for key, value in contentDict.items() if re.search('# cval', key)][0]
        minnhet = [value.strip() for key, value in contentDict.items() if re.search('min.nhet', key)][0]
        genome = [value.strip() for key, value in contentDict.items() if re.search('genome', key)][0]
        purity = [value.strip() for key, value in contentDict.items() if re.search('Purity', key)][0]
        ploidy = [value.strip() for key, value in contentDict.items() if re.search('Ploidy', key)][0]
        diplogr = [value.strip() for key, value in contentDict.items() if re.search('dipLogR', key)][0]
        loglik = [value.strip() for key, value in contentDict.items() if re.search('loglik', key)][0]
        flagsIdx = [i for i,x in enumerate(content) if 'flags' in x]
        flags = '|'.join(content[flagsIdx[0]+1:]).strip('#')
        output.write('\t'.join([sample, facets, snpnbhd, ndepth, puritycval, cval, minnhet, genome, purity, ploidy, diplogr, loglik, flags]) + '\n')

output.close()
