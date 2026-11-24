# Universality and collapse of completely positive factorizations

Verification scripts accompanying the paper by Minhao Liu (Texas Tech University). The scripts reconstruct the matrices and factorizations from the paper and compare all prescribed entries in the finite examples.

| Script | What it checks |
| --- | --- |
| `verify.py` | Four-bit moment recovery, the rational lift and projection of shortest factorizations, a four-cycle example in Shitov’s reduction, and the optional symbolic checks below. |
| `verify_paper_examples.py` | Principal deletions, parity and compatibility witnesses, and supplied rational factors. |
| `verify_rank_reductions.py` | Moment identities for tensors, orthogonal arrays, and Latin square completion examples. |

## Running the checks

The commands below use the Python standard library. From the repository directory, run:

```sh
python verify.py
python verify_paper_examples.py
python verify_rank_reductions.py
```

With NumPy and SciPy installed, the following command also checks Latin square examples by mixed-integer optimization.

```sh
python verify_rank_reductions.py --milp
```


With SymPy installed, run the arithmetic checks for the arithmetic examples in the manuscript. They verify irreducibility and the Galois group of the quintic, count its real roots, determine both fixed points of the fifth iterate of the Möbius map, and reconstruct both five-cycle Gram decompositions.

```sh
python verify.py --algebra
```

## Ordered private summands and Figure 1

Run `python check_private_order.py` to compare two shortest decompositions of the same rectangular matrix. It checks both complete lifts and the order of all 96 private summands using trace-normalized diagonal entries. Run `python make_figure_collapse.py` with Matplotlib to regenerate `fig_collapse.pdf`.

The optional Latin-square examples document related computations from earlier versions. The public code snapshot for this manuscript is [paper-2026-11-24](https://github.com/chineseschrodingercat/cp-factorization-verification/tree/paper-2026-11-24).
