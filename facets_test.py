import tempfile, subprocess, shutil, os, sys

TEST_TEMP_DIR = None
TEST_DATA_DIR = "/ifs/res/pwg/tests/cmo_testdata/facets/"


test_inputs = {"tumor_bam":"Chr22_hg19_TCGA-A8-A094-01A-11W-A019-09.tumor.bam",
        "chr22_vcf":"dbsnp_137.chr22.snp.positions.hg19.vcf.gz",
        "tumor_basecounts":"H_LS-A8-A094-01A-11W-A019-09-1.dat.gz",
        "normal_basecounts":"H_LS-A8-A094-10A-01W-A021-09-1.dat.gz",
        'maf':"TCGA-A8-A094-01A-11W-A019-09.maf",
        'facets_rdata_pairing':"facets_files.txt",
        'merged_counts':"countsMerged____TCGA-A8-A094-01A-11W-A019-09____TCGA-A8-A094-10A-01W-A021-09.dat.gz",
        'cncf_file':'H_LS-A8-A094-01A-11W-A019-09-1__H_LS-A8-A094-10A-01W-A021-09-1.cncf.txt',
        'big_maf_rdata_pairing':"rdata_big_mafanno/facets_data_file",
        'big_maf':"rdata_big_mafanno/Proj_04525_J___SOMATIC.vep.maf"
        }

for key, value in test_inputs.items():
    test_inputs[key]=os.path.join(TEST_DATA_DIR, "inputs", value)

expected_outputs = {"tumor_basecounts":"H_LS-A8-A094-01A-11W-A019-09-1.dat.gz",
        "mergeTN":"countsMerged____H_LS-A8-A094-01A-11W-A019-09-1__H_LS-A8-A094-10A-01W-A021-09-1.dat.gz",
        'seeded.seg':"H_LS-A8-A094-01A-11W-A019-09-1__H_LS-A8-A094-10A-01W-A021-09-1.seg",
        'seeded.seg.pcval':'H_LS-A8-A094-01A-11W-A019-09-1__H_LS-A8-A094-10A-01W-A021-09-1_hisens.seg',
        'ann_maf':"TCGA-A8-A094-01A-11W-A019-09.ann.maf",
        'gene_level_calls':"facets_gene_level_calls.txt",
        'big_ann_maf':"Proj_04525_J___SOMATIC.vep.facets_anno.maf",
        'arm_level_calls_em':"facets_arm_level_calls_em.txt",
        'arm_level_calls_cncf':"facets_arm_level_calls_cncf.txt"
        }

for key, value in expected_outputs.items():
    expected_outputs[key]=os.path.join(TEST_DATA_DIR, "expected_outputs", value)


DEV_NULL = open("/dev/null", "w")

def setup_module():
    global TEST_TEMP_DIR
    TEST_TEMP_DIR = tempfile.mkdtemp();

FACETS_SCRIPT = os.path.join(os.path.dirname(os.path.realpath(__file__)), "facets")

def teardown_module():
    if TEST_TEMP_DIR is not None:
        #shutil.rmtree(TEST_TEMP_DIR)
        pass

#def test_getbasecounts():
#    input_bam_file = test_inputs['tumor_bam']
#    input_vcf_file = test_inputs['chr22_vcf']
#    output_file = os.path.join(TEST_TEMP_DIR, "test_getbascounts_output")
#    cmd = ["cmo_getbasecounts", "--bam", input_bam_file,
#            "--out", output_file, "--sort_output", "--compress_output",  "--filter_improper_pair", "0",
#            "--vcf", input_vcf_file ]
#    rv = subprocess.call(cmd)
#    assert rv==0, "cmo_getbascounts threw shell error (program failed)"
#    expected_output= expected_outputs['tumor_basecounts']    
#    diff_cmd = ['diff', expected_output, output_file+".gz"]
#    rv = subprocess.call(diff_cmd)
#    assert rv==0, "getbasecounts output does not match expected output"

def test_mergeTN():
    input_tumor = test_inputs['tumor_basecounts']
    input_normal = test_inputs['normal_basecounts']
    test_output = os.path.join(TEST_TEMP_DIR, "test_countsmerged.gz")
    cmd = [FACETS_SCRIPT, "mergeTN", "-t", input_tumor, "-n",input_normal, "-o", test_output]
    rv = subprocess.call(cmd)
    expected_output= expected_outputs['mergeTN']
    shutil.copy(expected_output, TEST_TEMP_DIR)
    temp_expected_output = os.path.join(TEST_TEMP_DIR, os.path.basename(expected_output))
    gunzip_cmd = ["gunzip", temp_expected_output]
    subprocess.call(gunzip_cmd)
    gunzip_cmd = [ "gunzip", test_output]
    subprocess.call(gunzip_cmd)
    temp_expected_output= temp_expected_output.replace(".gz", "")
    test_output = test_output.replace(".gz", "")
    diff_cmd = ["diff", test_output, temp_expected_output]
    rv = subprocess.call(diff_cmd)
    assert rv==0, "cmo_facets mergeTN output does not match expected output, diff exit code: %s" % str(rv)

def test_facets():
    output_dir = os.path.join(TEST_TEMP_DIR)
    merged_count_input = test_inputs['merged_counts']
    #FIXME this test is hardcoded for our env...do more better...
    facets_cmd = [FACETS_SCRIPT, "doFacets", "-r=/opt/common/CentOS_6-dev/facets_lib/0.3.9/",
            "--seed=1587443596", 
            "-f", merged_count_input,
            "-t", "H_LS-A8-A094-01A-11W-A019-09-1__H_LS-A8-A094-10A-01W-A021-09-1",
            "-D", TEST_TEMP_DIR,
            "-c 200"]
    rv = subprocess.call(facets_cmd)
    assert rv==0, "facets failed to run :("
    expected_seg_output = expected_outputs['seeded.seg']
    test_seg_output = os.path.join(TEST_TEMP_DIR, "H_LS-A8-A094-01A-11W-A019-09-1__H_LS-A8-A094-10A-01W-A021-09-1.seg")
    diff_cmd = ["diff", expected_seg_output, test_seg_output]
    rv = subprocess.call(diff_cmd)
    assert rv==0, "facets test seg output differs from expected output- diff exit code: %s" % str(rv)

def test_facets_with_pcval():
    output_dir = os.path.join(TEST_TEMP_DIR)
    merged_count_input = test_inputs['merged_counts']
    facets_cmd = [FACETS_SCRIPT, "doFacets", "-r=/opt/common/CentOS_6-dev/facets_lib/0.3.9/",
            "--seed=1587443596", 
            "-f", merged_count_input,
            "-t", "H_LS-A8-A094-01A-11W-A019-09-1__H_LS-A8-A094-10A-01W-A021-09-1",
            "-D", TEST_TEMP_DIR,
            "--purity_cval=300",
            "-c=200"]
    print >>sys.stderr ," ".join(facets_cmd)
    rv = subprocess.call(facets_cmd)
    assert rv==0, "facets failed to run :("
    expected_seg_output = expected_outputs['seeded.seg.pcval']
    test_seg_output = os.path.join(TEST_TEMP_DIR, "H_LS-A8-A094-01A-11W-A019-09-1__H_LS-A8-A094-10A-01W-A021-09-1_hisens.seg")
    diff_cmd = ["diff", expected_seg_output, test_seg_output]
    rv = subprocess.call(diff_cmd)
    assert rv==0, "facets test seg output differs from expected output!"

def test_facets_maf():
    output_dir = os.path.join(TEST_TEMP_DIR)
    input_maf = test_inputs['maf']
    input_rdata_pairing = test_inputs['facets_rdata_pairing']
    facets_cmd = [FACETS_SCRIPT,
            "mafAnno",
            "-m", input_maf,
            "-f", input_rdata_pairing,
            "-o", os.path.join(TEST_TEMP_DIR, "TCGA-A8-A094-01A-11W-A019-09.ann.maf")]
    print " ".join(facets_cmd)
    rv = subprocess.call(facets_cmd)
    assert rv==0, "facets failed to run :("
    expected_annmaf_output = expected_outputs['ann_maf']
    test_seg_output = os.path.join(TEST_TEMP_DIR, "TCGA-A8-A094-01A-11W-A019-09.ann.maf")
    diff_cmd = ["diff", expected_annmaf_output, test_seg_output]
    rv = subprocess.check_call(diff_cmd)
    assert rv==0, "facets test seg output differs from expected output, diff exit code: %s" % str(rv)

def test_facets_big_maf():
    output_dir = os.path.join(TEST_TEMP_DIR)
    input_maf = test_inputs['big_maf']
    input_rdata_pairing = test_inputs['big_maf_rdata_pairing']
    facets_cmd = [FACETS_SCRIPT,
            "mafAnno",
            "-m", input_maf,
            "-f", input_rdata_pairing,
            "-o", os.path.join(TEST_TEMP_DIR, "big.ann.maf")]
    print " ".join(facets_cmd)
    rv = subprocess.call(facets_cmd)
    assert rv==0, "facets failed to run :("
    expected_annmaf_output = expected_outputs['big_ann_maf']
    test_seg_output = os.path.join(TEST_TEMP_DIR, "big.ann.maf")
    diff_cmd = ["diff", expected_annmaf_output, test_seg_output]
    rv = subprocess.check_call(diff_cmd)
    assert rv==0, "facets test seg output differs from expected output, diff exit code: %s" % str(rv)



def test_facets_gene_call():
    output_dir = os.path.join(TEST_TEMP_DIR)
    cncf_input = test_inputs['cncf_file']
    test_seg_output = os.path.join(TEST_TEMP_DIR, "facets_gene_level_calls.txt")
    facets_cmd = [FACETS_SCRIPT,
            "geneLevel",
            "-f", cncf_input,
            "-o",
            test_seg_output]
    print " ".join(facets_cmd)
    rv = subprocess.call(facets_cmd)
    assert rv==0, "facets failed to run :(, exit code %s" % rv
    expected_seg_output = expected_outputs['gene_level_calls']
    diff_cmd = ["diff", expected_seg_output, test_seg_output]
    rv = subprocess.call(diff_cmd)
    assert rv==0, "facets test seg output differs from expected output, diff exit code: %s" % str(rv)


def test_facets_arm_call_em():
    output_dir = os.path.join(TEST_TEMP_DIR)
    cncf_input = test_inputs['cncf_file']
    test_seg_output = os.path.join(TEST_TEMP_DIR, "facets_arm_level_calls_em.txt")
    facets_cmd = [FACETS_SCRIPT,
            "armLevel",
            "-f", cncf_input,
            "-o",
            test_seg_output,
            "-m",
            "em"]
    print " ".join(facets_cmd)
    rv = subprocess.call(facets_cmd)
    assert rv==0, "facets failed to run :(, exit code %s" % rv
    expected_seg_output = expected_outputs['arm_level_calls_em']
    diff_cmd = ["diff", expected_seg_output, test_seg_output]
    rv = subprocess.call(diff_cmd)
    assert rv==0, "facets test seg output differs from expected output, diff exit code: %s" % str(rv)


def test_facets_arm_call_cncf():
    output_dir = os.path.join(TEST_TEMP_DIR)
    cncf_input = test_inputs['cncf_file']
    test_seg_output = os.path.join(TEST_TEMP_DIR, "facets_arm_level_calls_cncf.txt")
    facets_cmd = [FACETS_SCRIPT,
            "armLevel",
            "-f", cncf_input,
            "-o",
            test_seg_output,
            "-m",
            "cncf"]
    print " ".join(facets_cmd)
    rv = subprocess.call(facets_cmd)
    assert rv==0, "facets failed to run :(, exit code %s" % rv
    expected_seg_output = expected_outputs['arm_level_calls_cncf']
    diff_cmd = ["diff", expected_seg_output, test_seg_output]
    rv = subprocess.call(diff_cmd)
    assert rv==0, "facets test seg output differs from expected output, diff exit code: %s" % str(rv)
