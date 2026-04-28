"""Dataset-driven regression tests for exact and heuristic solver quality."""

from __future__ import annotations

import networkx as nx
import pytest

from scanwidth import DAG, edge_scanwidth, node_scanwidth


def _graph_from_edges(edges: list[tuple[object, object]]) -> nx.DiGraph:
    """Build a DAG from explicit edge list."""
    graph = nx.DiGraph()
    graph.add_edges_from(edges)
    return graph


def _build_multi_root_1() -> nx.DiGraph:
    """Build explicit multi-root DAG #1."""
    return nx.DiGraph([(1, 3), (2, 3), (3, 4)])


REAL_NETWORK_EDGES = {
    "n01": [
        ("r", "94864"), ("r", "b"), ("b", "c"), ("c", "94864"), ("94864", "d"),
        ("d", "DBVPG_1373"), ("d", "L_1374"), ("d", "YJM978"), ("c", "e"),
        ("e", "UWOPS03.461.4"), ("e", "UWOPS05.217.3"), ("e", "UWOPS05.227.2"),
        ("c", "95044"), ("95044", "Y12"), ("b", "f"), ("f", "g"), ("g", "95044"),
        ("g", "Y9"), ("f", "h"), ("h", "i"), ("i", "K11"), ("i", "j"),
        ("j", "YPS606"), ("j", "YPS128"), ("i", "k"), ("k", "NCYC110"),
        ("k", "DBVPG6044"),
    ],
    "n03": [
        ("1", "Pan"), ("1", "2"), ("2", "3"), ("3", "4"), ("4", "Ccde"),
        ("4", "DIV_type10"), ("4", "5"), ("5", "DIII_type5"), ("3", "6"),
        ("6", "DIII_type4"), ("2", "7"), ("7", "8"), ("8", "9"), ("9", "10"),
        ("10", "5"), ("10", "weak_Dtype40"), ("9", "11"),
        ("11", "weak_Dtype420DAR"), ("11", "weak_Dtype421"), ("8", "DFV"),
        ("8", "12"), ("12", "DTO"), ("12", "RHDPSI"), ("7", "13"), ("13", "14"),
        ("14", "DAU_0"), ("14", "DAU_1"), ("13", "15"), ("15", "DAU_6"),
        ("15", "DAU_2"), ("7", "16"), ("16", "17"), ("17", "DNU"),
        ("17", "DIV_type5"), ("17", "RHD"), ("16", "18"), ("18", "6"),
        ("18", "19"), ("19", "DIV_type2"), ("19", "DIV_type3"),
        ("19", "DIV_type4"), ("19", "DNT"), ("19", "DWI"), ("19", "RHD_Ce"),
        ("16", "RHD_ce"), ("16", "RHD_deletion"), ("7", "DSF"), ("7", "DDN"),
        ("7", "DWN"),
    ],
    "n05": [
        ("r", "w1"), ("w1", "b1"), ("w1", "b2"), ("b1", "Actinobacteria"),
        ("b1", "b4"), ("r", "w2"), ("w2", "w3"), ("w2", "y1"), ("y1", "b2"),
        ("b2", "g1"), ("g1", "Firmicutes"), ("g1", "b4"), ("b4", "g3"),
        ("g3", "Proteobacteria"), ("g3", "g4"), ("g4", "Cyanobacteria"),
        ("g4", "g5"), ("g5", "p1"), ("p1", "p2"), ("p2", "Plants"),
        ("p2", "Animals"), ("y1", "y2"), ("y2", "Halobacteria"),
        ("w3", "Methanogens"), ("w3", "o1"), ("o1", "y2"), ("o1", "r1"),
        ("r1", "Eocytes"), ("r1", "p1"),
    ],
    "n07": [
        ("r", "AncestralModernHuman"), ("AncestralModernHuman", "African"),
        ("AncestralModernHuman", "BasalNonAfrican"), ("BasalNonAfrican", "h1"),
        ("h1", "v1"), ("v1", "v2"), ("v2", "European"), ("v2", "h2"),
        ("h2", "h3"), ("h3", "v3"), ("v3", "EE_NA"), ("v3", "h4"),
        ("h4", "Solomons"), ("v1", "h5"), ("h5", "v4"), ("v4", "h3"),
        ("v4", "AncestralOceanian"), ("AncestralOceanian", "v5"),
        ("v5", "Australian"), ("v5", "v6"), ("v6", "h4"), ("v6", "NewGuinea"),
        ("r", "v7"), ("v7", "v8"), ("v8", "Neanderthal"), ("v8", "v9"),
        ("v9", "h1"), ("v9", "h2"), ("v7", "v10"), ("v10", "h5"),
        ("v10", "Denisovan"),
    ],
    "n08": [
        ("r", "h1"), ("h1", "v1"), ("v1", "v2"), ("v2", "col"), ("v2", "gam"),
        ("v1", "h2"), ("h2", "ara"), ("r", "v3"), ("v3", "v4"), ("v4", "v5"),
        ("v5", "v6"), ("v6", "h1"), ("v6", "h2"), ("v5", "h3"), ("h3", "qua"),
        ("v4", "mel"), ("v7", "h3"), ("v3", "v7"), ("v7", "mer"),
    ],
    "n09": [
        ("r", "v1"), ("v1", "v2"), ("v2", "v3"), ("v3", "Gtricolor"),
        ("v3", "Gangelensis"), ("Gangelensis", "h1"), ("h1", "SantaAnaCanyon"),
        ("v2", "h2"), ("h2", "Gachilleafolia"), ("Gachilleafolia", "h3"),
        ("h3", "Gclivorum"), ("v1", "Gmillefoliata"), ("Gmillefoliata", "h3"),
        ("Gmillefoliata", "h4"), ("h4", "CapeMendocino"), ("r", "v4"),
        ("v4", "v5"), ("v5", "h2"), ("v5", "h5"), ("h5", "Gcapitata"),
        ("Gcapitata", "h4"), ("Gcapitata", "h1"), ("v4", "h5"),
    ],
    "n10": [
        ("r", "h1"), ("h1", "v1"), ("v1", "col"), ("v1", "h2"), ("h2", "gam"),
        ("r", "v2"), ("v2", "v3"), ("v3", "v4"), ("v4", "v5"), ("v5", "h1"),
        ("v5", "ara"), ("v4", "v6"), ("v6", "v7"), ("v7", "h2"), ("v7", "qua"),
        ("v6", "h3"), ("h3", "mer"), ("v3", "v8"), ("v8", "mel"), ("v8", "h4"),
        ("h4", "h3"), ("v2", "h4"),
    ],
    "n12": [
        ("UniversalDTU", "DTU1"), ("DTU1", "DTUI"), ("DTU1", "Hybrid12b"),
        ("Hybrid12b", "v1"), ("v1", "DTU2a"), ("DTU2a", "DTUIIa"),
        ("v1", "DTU2c"), ("DTU2c", "DTUIIc"), ("DTU2c", "Hybrid2b2c"),
        ("Hybrid2b2c", "v2"), ("v2", "DTU2d"), ("v2", "DTU2e"),
        ("UniversalDTU", "DTU2b"), ("DTU2b", "DTU2bp"), ("DTU2b", "Hybrid12b"),
        ("DTU2bp", "Hybrid2b2c"), ("DTU2bp", "DTUIIb"),
    ],
    "n14": [
        ("root", "Mbuti"), ("root", "NonAfrican"), ("NonAfrican", "v1"),
        ("v1", "EasternNonAfrican"), ("EasternNonAfrican", "h1"),
        ("h1", "Karitiana"), ("EasternNonAfrican", "Onge"), ("v1", "v2"),
        ("v2", "AncientNorthEurasian"), ("AncientNorthEurasian", "h1"),
        ("AncientNorthEurasian", "ANE"), ("ANE", "h2"), ("h2", "h3"),
        ("h3", "European"), ("ANE", "MA1"), ("v2", "WestEurasian"),
        ("WestEurasian", "WHG"), ("WHG", "Loschbour"), ("WHG", "h2"),
        ("WestEurasian", "EEF"), ("EEF", "h3"), ("EEF", "Stuttgart"),
        ("NonAfrican", "BasalEurasian"), ("BasalEurasian", "EEF"),
    ],
    "n15": [
        ("r", "Triticum"), ("r", "Aegilops"), ("Triticum", "h1"),
        ("Triticum", "v1"), ("Aegilops", "h1"), ("Aegilops", "v2"), ("h1", "v3"),
        ("v1", "Tmonococcum"), ("v1", "v4"), ("v2", "v5"), ("v2", "Aspeltoides"),
        ("v3", "Asharonensis"), ("v3", "v6"), ("v4", "Tuartu"), ("v4", "h2"),
        ("v5", "h2"), ("v5", "notIdentified"), ("v6", "Atauschii"),
        ("v6", "h3"), ("h2", "v7"), ("h3", "Taestivum"), ("v7", "h3"),
        ("v7", "Tturgidum"),
    ],
}


def _build_multi_root_2() -> nx.DiGraph:
    """Build explicit multi-root DAG #2."""
    return nx.DiGraph([(1, 4), (2, 4), (3, 5), (4, 6), (5, 6)])


def _build_multi_root_3() -> nx.DiGraph:
    """Build explicit multi-root DAG #3."""
    return nx.DiGraph([(1, 4), (2, 4), (2, 5), (3, 5), (4, 6), (5, 7)])


NETWORK_DATASET = [
    ("n01", _graph_from_edges(REAL_NETWORK_EDGES["n01"]), 2, 3),
    ("n03", _graph_from_edges(REAL_NETWORK_EDGES["n03"]), 2, 3),
    ("n05", _graph_from_edges(REAL_NETWORK_EDGES["n05"]), 3, 4),
    ("n07", _graph_from_edges(REAL_NETWORK_EDGES["n07"]), 3, 4),
    ("n08", _graph_from_edges(REAL_NETWORK_EDGES["n08"]), 3, 3),
    ("n09", _graph_from_edges(REAL_NETWORK_EDGES["n09"]), 4, 5),
    ("n10", _graph_from_edges(REAL_NETWORK_EDGES["n10"]), 4, 4),
    ("n12", _graph_from_edges(REAL_NETWORK_EDGES["n12"]), 2, 3),
    ("n14", _graph_from_edges(REAL_NETWORK_EDGES["n14"]), 3, 4),
    ("n15", _graph_from_edges(REAL_NETWORK_EDGES["n15"]), 3, 4),
    ("mr1", _build_multi_root_1(), 2, 2),
    ("mr2", _build_multi_root_2(), 2, 2),
    ("mr3", _build_multi_root_3(), 2, 2),
]

@pytest.mark.parametrize("name,graph,expected_nsw,expected_esw", NETWORK_DATASET)
def test_dataset_esw_solvers_against_baseline(
    name: str,
    graph: nx.DiGraph,
    expected_nsw: int,
    expected_esw: int,
) -> None:
    """All edge solvers are validated against the pinned ESW baseline."""
    dag = DAG(graph)
    exact_cases = [
        ("xp", {}),
        #("two_partition", {}),
        #("three_partition", {}),
    ]
    heuristic_cases = [
        ("greedy", {}),
        ("random", {"seed": 42}),
        ("cut_splitting", {}),
        ("simulated_annealing", {"max_iter": 100, "verbose": False, "seed": 42}),
    ]

    for algorithm, kwargs in exact_cases:
        value, _ = edge_scanwidth(dag, algorithm=algorithm, reduce=True, **kwargs)
        assert value == expected_esw, (
            f"{name}: {algorithm} should match exact esw={expected_esw}, got {value}"
        )

    for algorithm, kwargs in heuristic_cases:
        value, _ = edge_scanwidth(dag, algorithm=algorithm, reduce=True, **kwargs)
        assert value >= expected_esw, (
            f"{name}: heuristic {algorithm} returned {value} < exact esw={expected_esw}"
        )


@pytest.mark.parametrize("name,graph,expected_nsw,expected_esw", NETWORK_DATASET)
def test_dataset_nsw_solvers_against_baseline(
    name: str,
    graph: nx.DiGraph,
    expected_nsw: int,
    expected_esw: int,
) -> None:
    """All node solvers are validated against the pinned NSW baseline."""
    dag = DAG(graph)
    exact_cases = [
        ("xp", {}),
    ]
    heuristic_cases = [
        ("greedy", {}),
        ("random", {"seed": 42}),
        ("simulated_annealing", {"max_iter": 100, "verbose": False, "seed": 42}),
    ]

    for algorithm, kwargs in exact_cases:
        value, _ = node_scanwidth(dag, algorithm=algorithm, reduce=True, **kwargs)
        assert value == expected_nsw, (
            f"{name}: {algorithm} should match exact nsw={expected_nsw}, got {value}"
        )

    for algorithm, kwargs in heuristic_cases:
        value, _ = node_scanwidth(dag, algorithm=algorithm, reduce=True, **kwargs)
        assert value >= expected_nsw, (
            f"{name}: heuristic {algorithm} returned {value} < exact nsw={expected_nsw}"
        )
