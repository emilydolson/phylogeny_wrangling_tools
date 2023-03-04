import pandas as pd
import ALifeStdDev.phylogeny as phylodev
from pytest import approx
from .make_links import *
from .abstract_phylogenies import *


def verify_link_df(df):
    new_links = make_links(df)
    assert set(new_links.columns) == set(df["host"])
    assert set(new_links.index) == set(df["symbiont"])

    for h in new_links.columns:
        for s in new_links.index:

            sym_row = df[(df["host"] == h) & (df["symbiont"] == s)]
            if sym_row.empty:
                assert new_links.loc[s, h] == 0
            else:
                assert new_links.loc[s, h] == sym_row.iloc[0].at["count"]


def test_make_links():
    df = pd.read_csv("test_data/example_interaction_snapshot.csv")
    verify_link_df(df)
    # df = pd.read_csv("test_data/complex_interaction_snapshot.csv")
    # verify_link_df(df)


def test_abstract_phylogeny():
    original_df = phylodev.load_phylogeny_to_pandas_df("test_data/simple_phylogeny.csv")
    conversion_dict, df = processs_phylo(original_df, 20)

    # Check conversion dict
    for i in range(5):
        assert conversion_dict[i] == 0
    for i in range(6, 10):
        assert conversion_dict[i] == 0

    assert conversion_dict[10] == conversion_dict[11]
    assert conversion_dict[11] != 0
    assert conversion_dict[11] != conversion_dict[5]
    assert conversion_dict[5] != 0
    assert conversion_dict[12] != 0
    assert conversion_dict[12] != conversion_dict[5]
    assert conversion_dict[12] != conversion_dict[11]
    assert conversion_dict[13] != 0
    assert conversion_dict[13] != conversion_dict[5]
    assert conversion_dict[13] != conversion_dict[11]
    assert conversion_dict[13] != conversion_dict[12]
    assert conversion_dict[13] == conversion_dict[14]
    assert conversion_dict[13] == conversion_dict[15]

    # Check correspondence
    assert set(df["id"]) == set(conversion_dict.values())
    assert set(original_df.index) == set(conversion_dict.keys())

    # Check topology
    # Should have 3 rows
    assert (df.size / len(df.columns)) == 5
    assert df.loc[0, "ancestor_list"] == [None]
    assert df.loc[1, "ancestor_list"] == [0]
    assert df.loc[2, "ancestor_list"] == [0]
    assert df.loc[3, "ancestor_list"] == [3]
    assert df.loc[4, "ancestor_list"] == [3]
    assert df.loc[0, "edge_length"] == approx(0)
    assert df.loc[1, "edge_length"] == approx(.07)
    assert df.loc[2, "edge_length"] == approx(.07)

def test_abstract_phylogeny_complex():
    original_df = phylodev.load_phylogeny_to_pandas_df("test_data/complex_phylogeny.csv")
    conversion_dict, df = processs_phylo(original_df, 20)

    # Check correspondence
    assert set(df["id"]) == set(conversion_dict.values())
    assert set(original_df.index) == set(conversion_dict.keys())

    original_df = phylodev.load_phylogeny_to_pandas_df("test_data/complex_phylogeny_2.csv")
    conversion_dict, df = processs_phylo(original_df, 20)

    # Check correspondence
    assert set(df["id"]) == set(conversion_dict.values())
    assert set(original_df.index) == set(conversion_dict.keys())


def test_enrich():
    df = pd.read_csv("test_data/example_interaction_snapshot.csv")
    phylo_df = phylodev.load_phylogeny_to_pandas_df("test_data/simple_phylogeny.csv")
    df.columns = df.columns.str.replace(' ', '')

    df = enrich_interaction_df(df, phylo_df, phylo_df)

    assert set(df["host"]) == set(phylo_df.index)
    assert set(df["symbiont"]) == set(phylo_df.index)


def test_remove_excess_syms():
    df = pd.read_csv("test_data/example_interaction_snapshot.csv")
    df.columns = df.columns.str.replace(' ', '')
    phylo_df = phylodev.load_phylogeny_to_pandas_df("test_data/simple_phylogeny.csv")
    phylo_df.drop([2], inplace=True)
    df = remove_excess_symbionts(df, phylo_df)
    assert set(df["symbiont"]) == set([1, 3])


def test_integration():
    interaction_df = pd.read_csv("test_data/complex_interaction_snapshot.csv")
    interaction_df.columns = interaction_df.columns.str.replace(' ', '')

    sym_original_df = phylodev.load_phylogeny_to_pandas_df("test_data/complex_phylogeny.csv")
    sym_conversion_dict, sym_df = processs_phylo(sym_original_df, 20)
    host_original_df = phylodev.load_phylogeny_to_pandas_df("test_data/complex_phylogeny_2.csv")
    host_conversion_dict, host_df = processs_phylo(host_original_df, 20)

    interaction_df = remove_excess_symbionts(interaction_df, sym_original_df)
    interaction_df = enrich_interaction_df(interaction_df, sym_original_df, host_original_df)

    assert set(interaction_df["host"]) == set(host_original_df.index)
    assert set(interaction_df["symbiont"]) == set(sym_original_df.index)

    # print(interaction_df, host_df, host_conversion_dict)

    result = convert_interaction_labels(interaction_df, host_conversion_dict, sym_conversion_dict)

    assert set(result["host"]) == set(host_df["id"])
    assert set(result["symbiont"]) == set(sym_df["id"])

    new_links = make_links(result)
