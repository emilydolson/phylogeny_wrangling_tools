import pandas as pd
import ALifeStdDev.phylogeny as phylodev
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
    conversion_dict, df = processs_phylo("test_data/simple_phylogeny.csv", 20)

    # Check conversion dict
    for i in range(5):
        assert conversion_dict[i] == 0
    for i in range(6, 10):
        assert conversion_dict[i] == 0

    assert conversion_dict[10] == conversion_dict[11]
    assert conversion_dict[11] != 0
    assert conversion_dict[11] != conversion_dict[5]
    assert conversion_dict[5] != 0

    # Check correspondence
    assert set(df["id"]) == set(conversion_dict.values())
    assert set(original_df.index) == set(conversion_dict.keys())
