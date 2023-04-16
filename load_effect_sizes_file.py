import pandas as pd
def load_effect_sizes(effect_sizes):
    """Get effect size for each variant"""
    # Example: load_effect_sizes("resources/prs-scores-PGS002296.tsv")
    # effect_sizes = pd.read_csv(effect_sizes, sep = "\t")
    effect_sizes = effect_sizes[['chr', 'pos_hg38', 'risk_allele', 'weight']]
    effect_sizes['chr'] = effect_sizes['chr'].astype(str)
    effect_sizes['pos_hg38'] = effect_sizes['pos_hg38'].astype(int)
    effect_sizes['locus'] = effect_sizes['chr'].astype(str) + ':' + effect_sizes['pos_hg38'].astype(str) # pylint: disable=line-too-long
    effect_sizes_dict = effect_sizes[['locus', 'weight', 'risk_allele']].set_index('locus').T.to_dict() # pylint: disable=line-too-long

    return effect_sizes_dict
