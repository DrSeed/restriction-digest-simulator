# Self-contained demo: a plasmid digest on a virtual gel, and site frequencies vs expectation.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from digest import ENZYMES, find_sites, digest, expected_sites_per_mb

LADDER = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500]


def markov_dna(n, gc, rng, cpg_keep=1.0):
    # First-order Markov chain; cpg_keep < 1 depletes CpG like a mammalian genome.
    base = {"A": (1 - gc) / 2, "C": gc / 2, "G": gc / 2, "T": (1 - gc) / 2}
    letters = list("ACGT")
    seq = [rng.choice(letters, p=list(base.values()))]
    probs = {}
    for prev in letters:
        p = np.array([base[b] for b in letters])
        if prev == "C":
            p[2] *= cpg_keep
        probs[prev] = p / p.sum()
    draws = rng.random(n)
    for k in range(1, n):
        seq.append(letters[np.searchsorted(np.cumsum(probs[seq[-1]]), draws[k])])
    return "".join(seq)


def plasmid(rng):
    seq = list(markov_dna(5400, 0.5, rng))
    for site, pos in [("GAATTC", 400), ("GGATCC", 1450), ("AAGCTT", 2900), ("GAATTC", 4100)]:
        seq[pos:pos + len(site)] = list(site)
    return "".join(seq)


def mobility(bp):
    return -np.log10(bp)


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(5)
    pl = plasmid(rng)
    lanes = {"ladder": LADDER, "linearised (full length)": [len(pl)], "EcoRI": digest(pl, ["EcoRI"], True),
             "BamHI": digest(pl, ["BamHI"], True), "EcoRI+HindIII": digest(pl, ["EcoRI", "HindIII"], True),
             "EcoRI+BamHI+HindIII": digest(pl, ["EcoRI", "BamHI", "HindIII"], True)}

    models = {"random 50% GC": dict(gc=0.5), "random 40% GC": dict(gc=0.4),
              "40% GC, CpG depleted (mammal-like)": dict(gc=0.4, cpg_keep=0.22)}
    rows = []
    for label, kw in models.items():
        g = markov_dna(1_000_000, kw["gc"], rng, kw.get("cpg_keep", 1.0))
        for enz in ["EcoRI", "BamHI", "MspI", "NotI"]:
            rows.append({"sequence": label, "enzyme": enz, "site": ENZYMES[enz][0],
                         "observed_per_mb": len(find_sites(g, enz)),
                         "expected_per_mb_if_random": round(expected_sites_per_mb(ENZYMES[enz][0], kw["gc"]), 1)})
    table = pd.DataFrame(rows)
    table.to_csv("results/summary.csv", index=False)
    pd.DataFrame([{"lane": k, "fragments_bp": " ".join(map(str, v))} for k, v in lanes.items()]).to_csv("results/gel_lanes.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw={"width_ratios": [1, 1.4]})
    ax = axes[0]
    ax.set_facecolor("#1b1b1b")
    lo, hi = mobility(12000), mobility(300)
    for x, (name, frags) in enumerate(lanes.items()):
        for f in frags:
            y = (mobility(f) - lo) / (hi - lo)
            ax.add_patch(plt.Rectangle((x - 0.35, y - 0.006), 0.7, 0.012, color="#f5f5f5" if name != "ladder" else "#ffcc66"))
        if name == "ladder":
            for f in frags:
                ax.text(-0.5, (mobility(f) - lo) / (hi - lo), f"{f}", color="white", fontsize=7, ha="right", va="center")
    ax.set_xlim(-1.2, len(lanes) - 0.5)
    ax.set_ylim(1.02, -0.05)
    ax.set_xticks(range(len(lanes)), list(lanes), rotation=35, ha="right", fontsize=8)
    ax.set_yticks([])
    ax.set_title(f"Virtual 1% agarose gel, {len(pl)} bp plasmid")

    ax = axes[1]
    piv = table.pivot(index="enzyme", columns="sequence", values="observed_per_mb")[list(models)]
    piv.loc[["EcoRI", "BamHI", "MspI", "NotI"]].plot.bar(ax=ax, logy=True, rot=0)
    for i, enz in enumerate(["EcoRI", "BamHI", "MspI", "NotI"]):
        exp = expected_sites_per_mb(ENZYMES[enz][0], 0.4)
        ax.plot([i - 0.3, i + 0.3], [exp, exp], "k--", lw=1)
    ax.set_ylabel("sites per Mb (log scale)")
    ax.set_title("Observed sites; dashed = expected at 40% GC")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))
    print({k: v for k, v in lanes.items() if k != "ladder"})


if __name__ == "__main__":
    main()
