# Log

This file describes what has actually been done.

## COG0258 description

| Field                   | Information                                  |
| ----------------------- | -------------------------------------------- |
| **COG ID**              | **COG0258**                                  |
| **COG name**            | **5'-3' exonuclease Xni/ExoIX**              |
| **COG symbol**          | **ExoIX**                                    |
| **Functional category** | **L — Replication, recombination and repair**|
| **Gene/protein names**  | Xni; ExoIX / exonuclease IX                  |
| **PDB**                 | **3ZDA**                                     |
| **Literature**          | PubMed **19000038**, **23821668**            |

### Representative protein

The biologically characterized reference is **E. coli ExoIX**, encoded by
**xni**. The CDD page also provides representative protein sequences rather
than a single unique "characteristic protein" field. But there is an
experimentally characterized reference: Escherichia coli ExoIX/Xni.

## COG0258 analysis

Target `COG0258` was analysed with
[COGcollator](https://boabio.belozersky.msu.ru/COGcollator) with both filter off
and filter on.

![COGcollator-COG0258-filter-off](assets/01-COGcollator-COG0258-filter-off.png)
![COGcollator-COG0258-filter-on](assets/02-COGcollator-COG0258-filter-on.png)

I zoomed into the unfiltered `COG0258` graph in
[COGcollator](https://boabio.belozersky.msu.ru/COGcollator)

![03-COGcollator-filter-off-zoom-in](
    assets/03-COGcollator-filter-off-zoom-in.png
)

and put the selected sequences in
[DomainAnalyzer](https://boabio.belozersky.msu.ru/DomainAnalyser). The part with
non-`COG0258` proteins revealed the following picture:

![04-DomainAnalyzer-domain-overlap](assets/04-DomainAnalyzer-domain-overlap.png)
![05-DomainAnalyzer-tight-region](assets/05-DomainAnalyzer-tight-region.png)

This explains why `COG0749` was discarded almost entirely when filtering was
applied: it occurs frequently in fusion proteins as a second C-terminal domain.

The 7 remaining `COG0749` sequences also revealed no overlap:

![06-DomainAnalyzer-7-remaining](assets/06-DomainAnalyzer-7-remaining.png)

### Result

`COG0258` alone will define the initial COG-derived Xni/FEN homolog set.
No additional COG has been identified as a clearly related paralogous family
requiring inclusion.

## Relevant proteins extraction

The procedure aims to select

$$
\{\text{all COG0258 assignments}\}
\cap
\{\text{assignments from the 275 genomes}\}
$$

and then recover the corresponding amino-acid sequences from `COGorg24.faa`.

In order to accomplish that, I wrote a [parsing script](
    scripts/cog_parser/extract.py
) and ran it with

```sh
extract-cog \
    --cog COG0258 \
    --genomes data/27789-154191-1-SP.xlsx \
    --assignments data/large/cog-24.cog.csv \
    --sequences data/large/COGorg24.faa \
    --output data/COG0258 \
    --suppress-errors
```

which executed successfully

```sh
Representative genomes:       275
Selected COG assignments:     306
Genomes containing the COG:   273
Matched sequences:            306
Ambiguous FASTA records:      0
Duplicate FASTA matches:      0
Unmatched assignments:        0
Output directory:             ~/Desktop/Study/FBB/CourseWork2026/data/COG0258
```

and produced

```sh
data/COG0258/
├── COG0258_275_ambiguous.tsv
├── COG0258_275_raw.faa
├── COG0258_275_raw.tsv
└── COG0258_275_unmatched.tsv
```

- `COG0258_275_ambiguous.tsv` — FASTA records that could not be matched
  unambiguously;
- `COG0258_275_raw.faa` — corresponding amino-acid sequences;
- `COG0258_275_raw.tsv` — metadata for successfully matched assignments;
- `COG0258_275_unmatched.tsv` — selected COG assignments for which no
  sequence was recovered.

There are three source files, and each contributes different information:

```txt
27789-154191-1-SP.xlsx
        │
        │ Which genomes do we want?
        ▼
275 representative
    assemblies
        │
        │
        ├──────────────────────┐
        │                      │
        ▼                      │
  cog-24.cog.csv               │
        │                      │
        │ Which proteins       │
        │ belong to COG0258?   │
        ▼                      │
COG0258 assignments from ◄─────┘
    the 275 genomes
        │
        │ Which sequences correspond
        │ to those assignments?
        ▼
   COGorg24.faa
        │
        ▼
COG0258_275_raw.tsv
COG0258_275_raw.faa
```

### Result

The 2 files with filtered proteins sequences and assignments metadata

```sh
data/COG0258/
├── COG0258_275_raw.faa
└── COG0258_275_raw.tsv
```
