#!/usr/bin/env python


def get_protein_cluster_prevalence(df_input:pd.DataFrame):
    # Read Input
    genomes = sorted(df_input.iloc[:,0].unique())
    clusters = sorted(df_input.iloc[:,2].unique())

    # Create array
    A = np.zeros((len(genomes), len(clusters)), dtype=int)

    for _, (id_genome, id_protein, id_cluster) in df_input.iterrows():
        i = genomes.index(id_genome)
        j = clusters.index(id_cluster)
        A[i,j] += 1

    # Create output
    df_output = pd.DataFrame(A, index=genomes, columns=clusters)
    df_output.index.name = "id_genome"
    df_output.columns.name = "id_protein_cluster"

    return df_output