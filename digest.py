# Virtual restriction digests for linear or circular DNA.
import re

# name: (recognition site, cut offset on the top strand)
ENZYMES = {
    "EcoRI": ("GAATTC", 1), "BamHI": ("GGATCC", 1), "HindIII": ("AAGCTT", 1), "XhoI": ("CTCGAG", 1),
    "PstI": ("CTGCAG", 5), "EcoRV": ("GATATC", 3), "SmaI": ("CCCGGG", 3), "NdeI": ("CATATG", 2),
    "SpeI": ("ACTAGT", 1), "NotI": ("GCGGCCGC", 2), "MspI": ("CCGG", 1), "AluI": ("AGCT", 2),
}


def find_sites(seq, enzyme, circular=False):
    site, offset = ENZYMES[enzyme]
    seq = seq.upper()
    search = seq + seq[:len(site) - 1] if circular else seq
    # All the sites here are palindromic, so scanning the top strand finds every cut.
    return sorted({(m.start() + offset) % len(seq) for m in re.finditer(f"(?={site})", search)})


def digest(seq, enzymes, circular=False):
    """Fragment lengths after cutting with one or more enzymes."""
    cuts = sorted({c for e in enzymes for c in find_sites(seq, e, circular)})
    n = len(seq)
    if not cuts:
        return [n]
    if circular:
        return sorted(((cuts[(i + 1) % len(cuts)] - c) % n or n for i, c in enumerate(cuts)), reverse=True)
    edges = [0] + cuts + [n]
    return sorted((b - a for a, b in zip(edges, edges[1:]) if b > a), reverse=True)


def expected_sites_per_mb(site, gc=0.5):
    p = {"G": gc / 2, "C": gc / 2, "A": (1 - gc) / 2, "T": (1 - gc) / 2}
    prob = 1.0
    for b in site:
        prob *= p[b]
    return prob * 1e6
