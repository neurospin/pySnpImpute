##########################################################################
# NSAp - Copyright (C) CEA, 2016
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
##########################################################################

"""
TODO:
"""

from .info import __version__


# PACKAGE GLOBAL VARIABLES

# Chromosome names
ORDERED_AUTOSOMES = [str(c) for c in range(1, 23)]
X_PAR1, X_nonPAR, X_PAR2 = "X_PAR1", "X_nonPAR", "X_PAR2"
ORDERED_X_REGIONS = [X_PAR1, X_nonPAR, X_PAR2]
ORDERED_CHROMOSOMES = ORDERED_AUTOSOMES + ORDERED_X_REGIONS

# Sets
AUTOSOMES = set(ORDERED_AUTOSOMES)
X_REGIONS = set(ORDERED_X_REGIONS)
CHROMOSOMES = set(ORDERED_CHROMOSOMES)

# Common docs for argparse
REF_PANEL_DOC = (
    "Path to the table file describing the reference panel. "
    "It relates a chromosome name to 4 references files in Impute2 format:"
    " - reference haplotyples "
    " - reference legend "
    " - reference samples "
    " - recombination map "
    "The paths can be absolute or relative to the table file. "
    "File structure: one row per chromosome/X region, tab-separated fields: "
    "<CHROM>\t<REF HAP>\t<REF LEGEND>\t<REF SAMPLE>\t<RECOMBINATION MAP>"
    "Only the following chromosome names are accepted: "
    "'1', ..., '22', 'X_PAR1', 'X_nonPAR' and 'X_PAR2'"
)


# Positions to help with splitting the X chromosome in 3 regions:
# PAR1, nonPAR, PAR2
# The first key is the build version.
POSITIONS_OF_X_REGION_OF_BUILD = {
    "hg18": {"X_PAR1":   (1,           2709520),
             "X_nonPAR": (2709520 + 1, 154584238 - 1),
             "X_PAR2":   (154584238,   154913754)},
    "hg19": {"X_PAR1":   (60001,       2699520),
             "X_nonPAR": (2699520 + 1, 154931044 - 1),
             "X_PAR2":   (154931044,   155260560)},
    "hg38": {"X_PAR1":   (10001,       2781479),
             "X_nonPAR": (2781479 + 1, 155701383 - 1),
             "X_PAR2":   (155701383,   156030895)}
}

GENOMIC_BUILDS = set(POSITIONS_OF_X_REGION_OF_BUILD)


# Start and end positions of the centromere of each chromosome for each build
# Info extracted from the following UCSC files:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/cytoBand.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz
# http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
CENTROMERE_POSITIONS_OF_CHROMOSOME_OF_BUILD = {
    'hg18': {
        '1':       (124300000, 128000000),
        '2':        (93300000,  95700000),
        '3':        (91700000,  93200000),
        '4':        (50700000,  52400000),
        '5':        (47700000,  50500000),
        '6':        (60500000,  63400000),
        '7':        (59100000,  61100000),
        '8':        (45200000,  48100000),
        '9':        (51800000,  60300000),
        '10':       (40300000,  42100000),
        '11':       (52900000,  56400000),
        '12':       (35400000,  36500000),
        '13':       (16000000,  18400000),
        '14':       (15600000,  19100000),
        '15':       (17000000,  18400000),
        '16':       (38200000,  40700000),
        '17':       (22200000,  23200000),
        '18':       (16100000,  17300000),
        '19':       (28500000,  30200000),
        '20':       (27100000,  28400000),
        '21':       (12300000,  13200000),
        '22':       (11800000,  16300000),
        'X_nonPAR': (59500000,  65000000),
        'Y':        (11300000,  12500000)
    },
    'hg19': {
        '1':       (125000000, 128900000),
        '2':        (93300000,  96800000),
        '3':        (91000000,  93900000),
        '4':        (50400000,  52700000),
        '5':        (48400000,  50700000),
        '6':        (61000000,  63300000),
        '7':        (59900000,  61700000),
        '8':        (45600000,  48100000),
        '9':        (49000000,  50700000),
        '10':       (40200000,  42300000),
        '11':       (53700000,  55700000),
        '12':       (35800000,  38200000),
        '13':       (17900000,  19500000),
        '14':       (17600000,  19100000),
        '15':       (19000000,  20700000),
        '16':       (36600000,  38600000),
        '17':       (24000000,  25800000),
        '18':       (17200000,  19000000),
        '19':       (26500000,  28600000),
        '20':       (27500000,  29400000),
        '21':       (13200000,  14300000),
        '22':       (14700000,  17900000),
        'X_nonPAR': (60600000,  63000000),
        'Y':        (12500000,  13400000)
    },
    'hg38': {
        '1':       (123400000, 125100000),
        '2':        (93900000,  96000000),
        '3':        (90900000,  94000000),
        '4':        (50000000,  51800000),
        '5':        (48800000,  51400000),
        '6':        (59800000,  62600000),
        '7':        (60100000,  62100000),
        '8':        (45200000,  47200000),
        '9':        (43000000,  45500000),
        '10':       (39800000,  41600000),
        '11':       (53400000,  55800000),
        '12':       (35500000,  37800000),
        '13':       (17700000,  18900000),
        '14':       (17200000,  18200000),
        '15':       (19000000,  20500000),
        '16':       (36800000,  38400000),
        '17':       (25100000,  27400000),
        '18':       (18500000,  21500000),
        '19':       (26200000,  28100000),
        '20':       (28100000,  30400000),
        '21':       (12000000,  13000000),
        '22':       (15000000,  17400000),
        'X_nonPAR': (61000000,  63800000),
        'Y':        (10400000,  10600000)
    }
}


# Chromosome size for each build
# Info extracted from the following UCSC files
# ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/chromInfo.txt.gz
# ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz
# ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/chromInfo.txt.gz
SIZE_OF_CHROMOSOME_OF_BUILD = {
    'hg18': {
        '1':  247249719,
        '2':  242951149,
        '3':  199501827,
        '4':  191273063,
        '5':  180857866,
        '6':  170899992,
        '7':  158821424,
        '8':  146274826,
        '9':  140273252,
        '10': 135374737,
        '11': 134452384,
        '12': 132349534,
        '13': 114142980,
        '14': 106368585,
        '15': 100338915,
        '16':  88827254,
        '17':  78774742,
        '18':  76117153,
        '19':  63811651,
        '20':  62435964,
        '21':  46944323,
        '22':  49691432,
        'X':  154913754,
        'Y':   57772954
    },
    'hg19': {
        '1':  249250621,
        '2':  243199373,
        '3':  198022430,
        '4':  191154276,
        '5':  180915260,
        '6':  171115067,
        '7':  159138663,
        '8':  146364022,
        '9':  141213431,
        '10': 135534747,
        '11': 135006516,
        '12': 133851895,
        '13': 115169878,
        '14': 107349540,
        '15': 102531392,
        '16':  90354753,
        '17':  81195210,
        '18':  78077248,
        '19':  59128983,
        '20':  63025520,
        '21':  48129895,
        '22':  51304566,
        'X':  155270560,
        'Y':   59373566
    },
    'hg38': {
        '1':  248956422,
        '2':  242193529,
        '3':  198295559,
        '4':  190214555,
        '5':  181538259,
        '6':  170805979,
        '7':  159345973,
        '8':  145138636,
        '9':  138394717,
        '10': 133797422,
        '11': 135086622,
        '12': 133275309,
        '13': 114364328,
        '14': 107043718,
        '15': 101991189,
        '16':  90338345,
        '17':  83257441,
        '18':  80373285,
        '19':  58617616,
        '20':  64444167,
        '21':  46709983,
        '22':  50818468,
        'X':  156040895,
        'Y':   57227415
    }
}
