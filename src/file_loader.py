"""
Utilities to load graphs from disk into NetworkX with node positions.
Supports GraphML, GEXF, JSON (node-link), and EDGELIST/adjlist with optional node position attributes.

The loader tries to preserve or infer a `pos` attribute for each node as a tuple (x, y).
If positions are not present, a spring layout is computed and stored.
"""

import networkx as nx
import json
import os


def _try_parse_float(v):
    try:
        return float(v)
    except Exception:
        return None


def _ensure_pos_attr(g: nx.Graph):
    """Ensure each node has a `pos` attribute (tuple). If not present, compute spring layout."""
    has_pos = True
    for n, data in g.nodes(data=True):
        pos = data.get("pos")
        if pos is None:
            # try common keys
            x = data.get("x") or data.get("X") or data.get("cx") or data.get("px")
            y = data.get("y") or data.get("Y") or data.get("cy") or data.get("py")

            if x is not None and y is not None:
                xf = _try_parse_float(x)
                yf = _try_parse_float(y)
                if xf is not None and yf is not None:
                    g.nodes[n]["pos"] = (xf, yf)
                    continue

            # try a string like "x,y" or "(x, y)"
            for key in ["position", "pos", "coordinates"]:
                if key in data:
                    val = data[key]
                    if isinstance(val, str) and ("," in val or " " in val):
                        parts = val.replace("(", "").replace(")", "").replace("\n", " ").replace("\t", " ").split()
                        if len(parts) == 1 and "," in parts[0]:
                            parts = parts[0].split(",")
                        if len(parts) >= 2:
                            xf = _try_parse_float(parts[0].strip().strip(','))
                            yf = _try_parse_float(parts[1].strip().strip(','))
                            if xf is not None and yf is not None:
                                g.nodes[n]["pos"] = (xf, yf)
                                break
            else:
                has_pos = False

    if not has_pos:
        # compute layout
        pos = nx.spring_layout(g, seed=42)
        for n, p in pos.items():
            g.nodes[n]["pos"] = (float(p[0]), float(p[1]))

    return g


def load_graph(path: str) -> nx.Graph:
    """Load a graph from file and ensure `pos` on nodes.

    Supported formats by extension:
    - .graphml (recommended for editors like yEd)
    - .gexf
    - .json (node-link format)
    - .gpickle
    - .edgelist / .txt
    - .adjlist

    Returns a NetworkX Graph with node positions in attribute `pos`.
    """
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    ext = os.path.splitext(path)[1].lower()

    if ext == ".graphml":
        g = nx.read_graphml(path)
        # networkx may return a Graph with nodes as strings; try to convert numeric node ids if possible
        try:
            # attempt to cast node keys to int when possible
            mapping = {}
            for n in list(g.nodes()):
                try:
                    ni = int(n)
                    mapping[n] = ni
                except Exception:
                    pass
            if mapping:
                g = nx.relabel_nodes(g, mapping)
        except Exception:
            pass
        return _ensure_pos_attr(g)

    if ext == ".gexf":
        g = nx.read_gexf(path)
        return _ensure_pos_attr(g)

    if ext == ".json":
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        try:
            g = nx.node_link_graph(data)
        except Exception:
            # fallback to reading bare graph dict
            g = nx.node_link_graph(data)
        return _ensure_pos_attr(g)

    if ext in [".edgelist", ".txt"]:
        g = nx.read_edgelist(path)
        return _ensure_pos_attr(g)

    if ext == ".adjlist":
        g = nx.read_adjlist(path)
        return _ensure_pos_attr(g)

    # Unknown extension: try to infer via json
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        g = nx.node_link_graph(data)
        return _ensure_pos_attr(g)
    except Exception:
        pass

    raise ValueError(f"Unsupported file extension: {ext}")

