def aggregate_value_per_allele(df):
    df = df.copy()
    df["gene_id_collapsed"] = [
        ("_").join(gene_id.split("_")[:2]) for gene_id in df.gene_id
    ]
    df = df.drop(columns="gene_id").rename(columns={"gene_id_collapsed": "gene_id"})
    df = df.groupby("gene_id").agg("sum").reset_index()
    return df
