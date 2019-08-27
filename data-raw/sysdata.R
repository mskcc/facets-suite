suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(stringr)
    library(readr)
    library(purrr)
    library(usethis)
})

# Genome builds ---------------------------------------------------------------------------------------------------
# The code use to download the genome builds is commented, the resulting data is in the form of tribbles below.

# extract_genome_info = function(url) {
#     
#     # url can be a path to a single file or a list of urls if split by chromosome
#     build = map_dfr(url, function(x) {
#         read_tsv(x, col_names = c('bin', 'chrom', 'start', 'end', 'ix', 'n', 'size', 'type', 'bridge'))
#         }) %>% 
#         group_by(chrom) %>% 
#         summarize(
#             size = max(end),
#             centstart = min(start[type == 'centromere']),
#             centend = max(end[type == 'centromere']),
#             centromere = centstart + ((centend-centstart)/2)) %>% 
#         mutate(chrom = str_replace(chrom, 'chr', ''),
#                chrom = str_replace_all(chrom, c('X' = '23', 'Y' = '24')),
#                chrom = as.numeric(chrom)) %>% 
#         filter(chrom %in% seq(1, 24)) %>% 
#         arrange(chrom)
# }

# hg18 // pulled from UCSC
# hg18 = map(c(seq(1:22), 'X', 'Y'), ~paste0('http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chr', ., '_gap.txt.gz')) %>%
#     extract_genome_info()

hg18 = tibble::tribble(
    ~chrom,     ~size, ~centstart,  ~centend, ~centromere,
    1, 247249719,  121236957, 123476957,   122356957,
    2, 242951149,   91689898,  94689898,    93189898,
    3, 199501827,   90587544,  93487544,    92037544,
    4, 191273063,   49354874,  52354874,    50854874,
    5, 180857866,   46441398,  49441398,    47941398,
    6, 170899992,   58938125,  61938125,    60438125,
    7, 154817899,   58058273,  61058273,    59558273,
    8, 145403396,   43958052,  46958052,    45458052,
    9, 138336818,   47107499,  50107499,    48607499,
    10, 133577517,   39244941,  41624941,    40434941,
    11,  95942794,   51450781,  54450781,    52950781,
    12, 132349534,   34747961,  36142961,    35445461,
    13, 114142980,    1.6e+07,  17868000,    16934000,
    14, 106368585,   15070000,  18070000,    16570000,
    15,  96367672,   15260000,  18260000,    16760000,
    16,  88827254,   35143302,  36943302,    36043302,
    17,  78774742,   22187133,  22287133,    22237133,
    18,  73872808,   15400898,  16764896,    16082897,
    19,  63811651,   26923622,  29923622,    28423622,
    20,  60733814,   26267569,  28033230,  27150399.5,
    21,  43507092,   10260000,  13260000,    11760000,
    22,  49691432,   11330000,  14330000,    12830000,
    23, 148832720,   58598737,  61598737,    60098737,
    24,  57377044,   11253954,  12308578,    11781266
)

# hg19 // pulled from UCSC
# hg19 = extract_genome_info('http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz')
hg19 = tibble::tribble(
    ~chrom,     ~size, ~centstart,  ~centend, ~centromere,
    1, 249250621,  121535434, 124535434,   121535434,
    2, 243199373,   92326171,  95326171,    92326171,
    3, 198022430,   90504854,  93504854,    90504854,
    4, 191154276,   49660117,  52660117,    49660117,
    5, 180915260,   46405641,  49405641,    46405641,
    6, 171115067,   58830166,  61830166,    58830166,
    7, 159138663,   58054331,  61054331,    58054331,
    8, 146364022,   43838887,  46838887,    43838887,
    9, 141213431,   47367679,  50367679,    47367679,
    10, 135534747,   39254935,  42254935,    39254935,
    11, 135006516,   51644205,  54644205,    51644205,
    12, 133851895,   34856694,  37856694,    34856694,
    13, 115169878,    1.6e+07,   1.9e+07,     1.6e+07,
    14, 107349540,    1.6e+07,   1.9e+07,     1.6e+07,
    15, 102531392,    1.7e+07,     2e+07,     1.7e+07,
    16,  90354753,   35335801,  38335801,    35335801,
    17,  79759049,   22263006,  25263006,    22263006,
    18,  78077248,   15460898,  18460898,    15460898,
    19,  59128983,   24681782,  27681782,    24681782,
    20,  63025520,   26369569,  29369569,    26369569,
    21,  48129895,   11288129,  14288129,    11288129,
    22,  51304566,    1.3e+07,   1.6e+07,     1.3e+07,
    23, 155270560,   58632012,  61632012,    58632012,
    24,  59373566,   10104553,  13104553,    10104553
)

# hg38 // pulled from UCSC
# hg38 = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz' # note: centromere info not longer in this type of file, but in centromeres.txt.gz
# hg38_centromeres = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz' 
# hg38_centromeres = read_tsv(grch38_centromeres, col_names = c('bin', 'chrom', 'start', 'end', 'name')) %>% 
#     group_by(chrom) %>% 
#     summarize(centstart = min(start),
#               centend = max(end),
#               centromere = centstart + ((centend-centstart)/2)) %>% 
#     mutate(chrom = str_replace(chrom, 'chr', ''),
#            chrom = str_replace_all(chrom, c('X' = '23', 'Y' = '24')),
#            chrom = as.numeric(chrom))
# 
# hg38 = extract_genome_info(hg38) %>% 
#     select(-matches('cent')) %>% 
#     left_join(., hg38_centromeres, by = 'chrom')

hg38 = tibble::tribble(
    ~chrom,     ~size, ~centstart,  ~centend, ~centromere,
    1, 248956422,  122026459, 124932724, 123479591.5,
    2, 242193529,   92188145,  94090557,    93139351,
    3, 198295559,   90772458,  93655574,    92214016,
    4, 190214555,   49712061,  51743951,    50728006,
    5, 181538259,   46485900,  50059807,  48272853.5,
    6, 170805979,   58553888,  59829934,    59191911,
    7, 159345973,   58169653,  61528020,  59848836.5,
    8, 145138636,   44033744,  45877265,  44955504.5,
    9, 138394717,   43389635,  45518558,  44454096.5,
    10, 133797422,   39686682,  41593521,  40640101.5,
    11, 135086622,   51078348,  54425074,    52751711,
    12, 133275309,   34769407,  37185252,  35977329.5,
    13, 114364328,    1.6e+07,  18051248,    17025624,
    14, 107043718,    1.6e+07,  18173523,  17086761.5,
    15, 101991189,   17083673,  19725254,  18404463.5,
    16,  90338345,   36311158,  38265669,  37288413.5,
    17,  83257441,   22813679,  26616164,  24714921.5,
    18,  80373285,   15460899,  20861206,  18161052.5,
    19,  58617616,   24498980,  27190874,    25844927,
    20,  64444167,   26436232,  30038348,    28237290,
    21,  46709983,   10864560,  12915808,    11890184,
    22,  50818468,   12954788,  15054318,    14004553,
    23, 156040895,   58605579,  62412542,  60509060.5,
    24,  57227415,   10316944,  10544039,  10430491.5
)

# Gene positions --------------------------------------------------------------------------------------------------
genes_hg19 = fread('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh37_mapping/gencode.v29lift37.basic.annotation.gtf.gz',
                   header = F, skip = 'chr1',
                   col.names = c('chrom', 'source', 'type', 'start', 'end', 'na_1', 'strand', 'na_2', 'info')) %>%
    filter(type == 'gene', info %like% 'protein_coding') %>%
    mutate(chrom = str_replace(chrom, 'chr', ''),
           gene = str_extract(info, '(?<=gene_name ")[A-Za-z0-9\\.\\-]+(?=";)')) %>%
    group_by(gene, chrom) %>%
    summarize(start = min(start),
              end = max(end)) %>%
    ungroup()

genes_hg38 = fread('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.basic.annotation.gtf.gz',
                   header = F, skip = 'chr1',
                   col.names = c('chrom', 'source', 'type', 'start', 'end', 'na_1', 'strand', 'na_2', 'info')) %>%
    filter(type == 'gene', info %like% 'protein_coding') %>%
    mutate(chrom = str_replace(chrom, 'chr', ''),
           gene = str_extract(info, '(?<=gene_name ")[A-Za-z0-9\\.\\-]+(?=";)')) %>%
    group_by(gene, chrom) %>%
    summarize(start = min(start),
              end = max(end)) %>%
    ungroup()

# use_data(hg18, hg19, hg38, genes_hg19, genes_hg38, internal = T, overwrite = T)

# Copy-number states ----------------------------------------------------------------------------------------------
copy_number_states = tibble::tribble(
    ~wgd, ~mcn, ~lcn, ~numeric_call,          ~call,
    # No genome doubling
    FALSE,    0,    0,          -2,              'HOMDEL',
    FALSE,    1,    0,          -1,             'HETLOSS',
    FALSE,    2,    0,          -1,               'CNLOH',
    FALSE,    3,    0,           1,        'CNLOH & GAIN',
    FALSE,    4,    0,           1,        'CNLOH & GAIN',
    FALSE,    5,    0,           2,           'AMP (LOH)',
    FALSE,    6,    0,           2,           'AMP (LOH)',
    FALSE,    1,    1,           0,             'DIPLOID',
    FALSE,    2,    1,           1,                'GAIN',
    FALSE,    3,    1,           1,                'GAIN',
    FALSE,    4,    1,           2,                 'AMP',
    FALSE,    5,    1,           2,                 'AMP',
    FALSE,    6,    1,           2,                 'AMP',
    FALSE,    2,    2,           1,          'TETRAPLOID',
    FALSE,    3,    2,           2,                 'AMP',
    FALSE,    4,    2,           2,                 'AMP',
    FALSE,    5,    2,           2,                 'AMP',
    FALSE,    6,    2,           2,                 'AMP',
    FALSE,    3,    3,           2,      'AMP (BALANCED)',
    FALSE,    4,    3,           2,                 'AMP',
    FALSE,    5,    3,           2,                 'AMP',
    FALSE,    6,    3,           2,                 'AMP',
    # With genome doubling
    TRUE,    0,    0,          -2,              'HOMDEL',
    TRUE,    1,    0,          -1, 'LOSS BEFORE & AFTER',
    TRUE,    2,    0,          -1,         'LOSS BEFORE',
    TRUE,    3,    0,          -1, 'CNLOH BEFORE & LOSS',
    TRUE,    4,    0,          -1,        'CNLOH BEFORE',
    TRUE,    5,    0,           1, 'CNLOH BEFORE & GAIN',
    TRUE,    6,    0,           2,           'AMP (LOH)',
    TRUE,    1,    1,          -1,   'DOUBLE LOSS AFTER',
    TRUE,    2,    1,          -1,          'LOSS AFTER',
    TRUE,    3,    1,          -1,         'CNLOH AFTER',
    TRUE,    4,    1,           1,         'LOSS & GAIN',
    TRUE,    5,    1,           2,                 'AMP',
    TRUE,    6,    1,           2,                 'AMP',
    TRUE,    2,    2,           0,          'TETRAPLOID',
    TRUE,    3,    2,           1,                'GAIN',
    TRUE,    4,    2,           2,                 'AMP',
    TRUE,    5,    2,           2,                 'AMP',
    TRUE,    6,    2,           2,                 'AMP',
    TRUE,    3,    3,           2,      'AMP (BALANCED)',
    TRUE,    4,    3,           2,                 'AMP',
    TRUE,    5,    3,           2,                 'AMP',
    TRUE,    6,    3,           2,                 'AMP'
) %>% mutate(map_string := paste(wgd, mcn, lcn, sep = ':'))

use_data(hg18, hg19, hg38, genes_hg19, genes_hg38, copy_number_states, internal = T, overwrite = T)
