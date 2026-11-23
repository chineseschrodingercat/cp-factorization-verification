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


With SymPy installed, run the arithmetic checks for Remark 3.3. They verify irreducibility and the Galois group of the quintic, count its real roots, determine both fixed points of the fifth iterate of the Möbius map, and reconstruct both five-cycle Gram decompositions.

```sh
python verify.py --algebra
```
