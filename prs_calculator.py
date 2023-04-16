"""Script for calculation PRS in a additive fashion (sum(beta * genotype))"""
import pandas as pd

def calc_vcf_prs(vcf_file: str, risk_alleles_weight: dict) -> dict:
    """
    Calculate PRS of each sample by multiplying the effect_size (weight) by the number of risk alleles
    
    Args:
        vcf_filename (str): path to the VCF file;
        risk_alleles_weight (dict): dictionary as follows

        {'chr1:959139': {'weight': 0.05, 'risk_allele': 'G'},
        'chr1:1127258': {'weight': -0.016, 'risk_allele': 'C'},
        'chr1:1748780': {'weight': 0.021, 'risk_allele': 'G'},
        'chr1:2115499': {'weight': -0.019, 'risk_allele': 'G'}} 
    
    Returns:
        dict: as follows
        {'Sample_01': 2.3,
        'Sample_02': -4.5,
        'Sample_03': -2.8,
        'Sample_04': 3.6,
        'Sample_05': 5.17
        }   
    """
    

    got_header = False
    prs_sum = []
    # with open(vcf_filename, 'r') as vcf_file:
    # with vcf_filename as vcf_file:
    for raw_line in vcf_file:
        line = raw_line.decode("utf-8").split('\t')

        # Skip commentaries
        if raw_line[0:2].decode("utf-8") == '##':
            continue

        if not got_header:
            # get headers
            vcf_names = [x.strip() for x in line]
            vcf_names[0] = 'CHROM'


            index_chr = vcf_names.index('CHROM')
            index_pos = vcf_names.index('POS')

            got_header = True

        else:
            prs_sum = process_line(
                line, vcf_names, risk_alleles_weight,
                index_chr, index_pos, prs_sum)


    sample_names = vcf_names[9::]
    return dict(zip(sample_names, [ round(x, 5) for x in prs_sum ]))

def process_line(line, vcf_names, risk_alleles_weight, index_chr, index_pos, prs_sum): # pylint: disable=line-too-long
    """Calculate PRS sum for all samples"""
    locus = line[index_chr]+":"+line[index_pos]

    try:
        risk_variant_locus = risk_alleles_weight[locus]['weight']
        risk_allele = risk_alleles_weight[locus]['risk_allele']

        # get risk_allele index
        risk_allele_index = get_effect_risk_index(line, risk_allele, vcf_names)

        n_samples = len(vcf_names)

        sample_effects = [
            sum_risk_allele_index(
                sample_genotype = x.split(":")[0],
                risk_allele_index = risk_allele_index,
                risk_variant_locus = risk_variant_locus) for x in line[9:n_samples]
        ]

        return add_sample_effect(sample_effects, prs_sum)

    except KeyError:
        #print(f"Variant {locus} present in VCF is not present in PRS file")

        return prs_sum

def add_sample_effect(values_to_add, prs_sum):
    """Create sum of PRS of each sample"""
    if prs_sum == []:
        new_sum = values_to_add
    else:
        new_sum = [x + y for x, y in zip(values_to_add, prs_sum)]

    return new_sum


def sum_risk_allele_index(sample_genotype, risk_allele_index, risk_variant_locus):
    """Calc additive model"""
    return sample_genotype.count(str(risk_allele_index)) * risk_variant_locus

def get_effect_risk_index(line, risk_allele, vcf_names):
    """Get GT index for risk allele"""
    try:
        risk_allele_index = line[vcf_names.index('ALT')].split(",").index(risk_allele) + 1
    except ValueError:
        risk_allele_index = -1

    # if risk allele is are the reference allele
    if line[vcf_names.index('REF')] == risk_allele:
        risk_allele_index =  0

    return risk_allele_index
