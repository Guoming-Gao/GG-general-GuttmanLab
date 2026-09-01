from spen_metagene.clap import build_index, assign


def test_smallest_overlapping_gene_is_selected():
    models={
      "large":{"chrom":"chr1","strand":"+","start":100,"end":1000},
      "small":{"chrom":"chr1","strand":"+","start":200,"end":300},
    }
    idx=build_index(models,0)
    assert assign(idx,"chr1","+",250)=="small"
    assert assign(idx,"chr1","-",250) is None

