yaomatchromosome_dict = {
    'chr1_hap1': 'chr1',
    'chr2_hap1': 'chr2',
    'chr3_hap1': 'chr3',
    'chr4_hap1': 'chr4',
    'chr5_hap1': 'chr5',
    'chr6_hap1': 'chr6',
    'chr7_hap1': 'chr7',
    'chr8_hap1': 'chr8',
    'chr9_hap1': 'chr9',
    'chr10_hap1': 'chr10',
    'chr11_hap1': 'chr11',
    'chr12_hap1': 'chr12',
    'chr13_hap1': 'chr13',
    'chr14_hap1': 'chr14',
    'chr15_hap1': 'chr15',
    'chr16_hap1': 'chr16',
    'chr17_hap1': 'chr17',
    'chr18_hap1': 'chr18',
    'chr19_hap1': 'chr19',
    'chr20_hap1': 'chr20',
    'chr21_hap1': 'chr21',
    'chr22_hap1': 'chr22',
    'chrX_hap1': 'chrX'
}

yaopatchromosome_dict = {
    'chr1_hap2': 'chr1',
    'chr2_hap2': 'chr2',
    'chr3_hap2': 'chr3',
    'chr4_hap2': 'chr4',
    'chr5_hap2': 'chr5',
    'chr6_hap2': 'chr6',
    'chr7_hap2': 'chr7',
    'chr8_hap2': 'chr8',
    'chr9_hap2': 'chr9',
    'chr10_hap2': 'chr10',
    'chr11_hap2': 'chr11',
    'chr12_hap2': 'chr12',
    'chr13_hap2': 'chr13',
    'chr14_hap2': 'chr14',
    'chr15_hap2': 'chr15',
    'chr16_hap2': 'chr16',
    'chr17_hap2': 'chr17',
    'chr18_hap2': 'chr18',
    'chr19_hap2': 'chr19',
    'chr20_hap2': 'chr20',
    'chr21_hap2': 'chr21',
    'chr22_hap2': 'chr22',
    'chrY_hap2': 'chrY'
}

hg38chromosome_dict = nc_dict = {
    'NC_000001.11': 'chr1',
    'NC_000002.12': 'chr2',
    'NC_000003.12': 'chr3',
    'NC_000004.12': 'chr4',
    'NC_000005.10': 'chr5',
    'NC_000006.12': 'chr6',
    'NC_000007.14': 'chr7',
    'NC_000008.11': 'chr8',
    'NC_000009.12': 'chr9',
    'NC_000010.11': 'chr10',
    'NC_000011.10': 'chr11',
    'NC_000012.12': 'chr12',
    'NC_000013.11': 'chr13',
    'NC_000014.9': 'chr14',
    'NC_000015.10': 'chr15',
    'NC_000016.10': 'chr16',
    'NC_000017.11': 'chr17',
    'NC_000018.10': 'chr18',
    'NC_000019.10': 'chr19',
    'NC_000020.11': 'chr20',
    'NC_000021.9': 'chr21',
    'NC_000022.11': 'chr22',
    'NC_000023.11': 'chrX',
    'NC_000024.10': 'chrY',
    'NC_012920.1': 'chrM'
}

# 还可以创建反向映射（如果需要）
yaomatchromosome_dict_rev = {v: k for k, v in yaomatchromosome_dict.items()}
yaopatchromosome_dict_rev = {v: k for k, v in yaopatchromosome_dict.items()}
hg38chromosome_dict_rev = {v: k for k, v in hg38chromosome_dict.items()}
