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
