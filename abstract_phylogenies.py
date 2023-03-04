# Process host, symbiont, and interaction snapshots into abstracted phylogenies
# and an updated interaction snapshot

import sys
import math
import pandas as pd
import ALifeStdDev.phylogeny as phylodev
import networkx as nx


def processs_phylo(df, num_bins=1000, outfilename="compressed.csv"):

    min_val = -1
    max_val = 1

    df["Bin"] = df["info"].apply(
        lambda x: math.floor((num_bins - 1) * (x - min_val)/max_val))

    conversion_dict = {}

    g = phylodev.pandas_df_to_networkx(df)

    abstract_g = phylodev.abstract_asexual_phylogeny(g, ["Bin"])

    for node in abstract_g.nodes:
        members = abstract_g.nodes[node]["members"]
        original = min(members, key=lambda k: members[k]["origin_time"])
        abstract_g.nodes[node]["origin"] = original
        for m in members:
            conversion_dict[m] = node
        # abstract_g.nodes[node]["original_info"] = members[original]["info"]

    for node in abstract_g.nodes:
        original = abstract_g.nodes[node]["origin"]
        parent = list(abstract_g.predecessors(node))

        if len(parent) == 0:
            abstract_g.nodes[node]["edge_length"] = 0
            continue

        parent_original = abstract_g.nodes[parent[0]]["origin"]
        path = nx.shortest_path(g, source=parent_original, target=original)
        dist = 0
        curr = g.nodes[parent_original]["info"]
        for n in path:
            dist += abs(curr - g.nodes[n]["info"])
            curr = g.nodes[n]["info"]
        abstract_g.nodes[node]["edge_length"] = dist

    final_df = phylodev.networkx_to_pandas_df(abstract_g, {"edge_length":"edge_length"})
    final_df["taxon_label"] = final_df["id"]
    final_df.to_csv(outfilename, index=False)
    return conversion_dict, final_df


def convert_interaction_labels(df, host_dict, sym_dict):
    df.loc[:, "host"] = df.loc[:, "host"].apply(lambda x: host_dict[x])
    df.loc[:, "symbiont"] = df.loc[:, "symbiont"].apply(lambda x: sym_dict[x])

    agg_functions = {"host": "first",
                     "symbiont": "first",
                     "count": "sum"
                     }

    return df.groupby(by=["host", "symbiont"]).aggregate(agg_functions)


def enrich_interaction_df(interaction_df, sym_df, host_df):
    missing_syms = set(sym_df.index) - set(interaction_df["symbiont"])
    missing_hosts = set(host_df.index) - set(interaction_df["host"])

    token_sym = sym_df.index[0]
    token_host = host_df.index[0]

    for sym in missing_syms:
        new_row = pd.Series({"host": token_host, "symbiont": sym, "count": 0})
        interaction_df = pd.concat([interaction_df, new_row.to_frame().T], ignore_index=True)

    for host in missing_hosts:
        new_row = pd.Series({"host": host, "symbiont": token_sym, "count": 0})
        interaction_df = pd.concat([interaction_df, new_row.to_frame().T], ignore_index=True)
  
    return interaction_df


def remove_excess_symbionts(interaction_df, sym_phylo):
    return interaction_df[interaction_df["symbiont"].apply(lambda x: x in sym_phylo.index)]


def main():
    if len(sys.argv) <= 3:
        print("Usage: python abstract_phylogenies.py [Symbiont Phylogeny Snapshot file] [Host Phylogeny Snapshot file] [Interaction Snapshot file] [resolution (optional)]")

    sym_filename = sys.argv[1]
    host_filename = sys.argv[2]
    interaction_filename = sys.argv[3]
    num_bins = 1000
    if len(sys.argv) > 4:
        num_bins = int(sys.argv[4])
    
    sym_df = phylodev.load_phylogeny_to_pandas_df(sym_filename)
    host_df = phylodev.load_phylogeny_to_pandas_df(host_filename)

    sym_conversion_dict = processs_phylo(sym_df, num_bins)
    host_conversion_dict = processs_phylo(host_df, num_bins)

    interaction_df = pd.read_csv(interaction_filename)
    interaction_df.columns = interaction_df.columns.str.replace(' ', '')
    interaction_df = enrich_interaction_df(interaction_df, sym_df, host_df)
    interaction_df = remove_excess_symbionts(interaction_df, sym_df)

    result = convert_interaction_labels(interaction_df, host_conversion_dict, sym_conversion_dict)
    result.to_csv('abstract_interactions.csv', index=False)


if __name__ == "__main__":
    main()