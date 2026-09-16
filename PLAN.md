# Research Plan

## Definition of the protein group under study

1. Use COG0258 as the initial group of homologous proteins.
2. Analyze COG0258 using COGcollator.
3. Identify closely related COGs and potential paralogous groups.
4. Determine which COG(s) should be included in the initial dataset.
5. Establish criteria for inclusion of proteins in the studied group.

## Collection of sequences from the COG database

1. Download the required COG database files:
    - cog-24.cog.csv;
    - cog-24.org.csv;
    - COGorg24.faa.
2. Extract proteins assigned to the selected COG(s).
Restrict the dataset to the set of 275 representative prokaryotic genomes
*(see "27789-154191-1-SP.xlsx", TableS2 to learn more)*.
3. Collect the following metadata for each protein:
    - protein ID;
    - locus tag;
    - organism name;
    - taxonomic ID;
    - genome assembly;
    - COG assignment;
    - COG-footprint coordinates;
    - COG assignment quality.
This follows the structure of the current COG database files described in the
supplied materials. In particular, cog-24.cog.csv provides COG assignments,
footprint coordinates, and assignment-quality information.

## Preparation and quality control of the sequence dataset

1. Assess the quality of the extracted sequences.
2. Check:
    - completeness of COG-footprints;
    - presence of fragmented proteins;
    - anomalous sequence lengths;
    - possible annotation errors;
    - partial or low-quality COG assignments.
3. Remove clearly incorrect or artifactual sequences.
4. Retain ambiguous cases separately for manual examination.
5. Generate a curated COG-derived sequence dataset.

## Identification of experimentally characterized homologs

1. Search for homologous proteins using:
    - BLAST against PDB;
    - BLAST against Swiss-Prot.
2. Identify proteins with:
    - experimentally determined structures;
    - relevant experimental functional data.
3. Verify their COG/domain assignment using DomainAnalyser.
4. Exclude experimentally characterized proteins that do not belong to the
studied homologous group.
5. Add confirmed homologs to the sequence dataset.

This is also the workflow recommended in the supplied material: search
PDB/Swiss-Prot for experimentally studied homologs and then verify their COG
assignment before incorporating them into the analysis.

## Identification of Zn-binding residues

1. For experimentally and structurally characterized homologs:
    1. Collect relevant PDB structures and literature data.
    2. Record:
        - PDB ID;
        - protein/chain identifier;
        - bound metal ions;
        - Zn²⁺ ions, where present;
        - residues involved in Zn²⁺ coordination;
        - catalytic residues;
        - other functionally relevant active-site residues.
2. Create a reference table containing:
    - protein;
    - organism;
    - PDB ID;
    - Zn-binding residues;
    - catalytic residues;
    - relevant experimental evidence.

## Multiple sequence alignment (MSA)

1. Combine:
    - curated COG-derived sequences;
    - confirmed experimentally characterized homologs.
2. Extract the homologous Xni/FEN regions where necessary.
3. Perform multiple sequence alignment (MSA).
4. Assess MSA quality.
5. Identify and remove:
    - poorly aligned sequences;
    - remaining sequence fragments;
    - unreliable alignment regions;
    - residual artifacts.
6. Generate the final filtered MSA.

## Mapping of Zn-binding positions onto the MSA

1. Map experimentally identified Zn-binding residues to the corresponding MSA
columns.
2. Map other relevant catalytic and active-site residues.
3. Annotate the MSA with:
    - Zn-binding positions;
    - catalytic positions;
    - other conserved active-site positions.
4. Determine the amino acid state of each Zn-binding position in every sequence.
5. Evaluate the conservation of individual Zn-binding residues and of the
Zn-binding site as a whole.

## Phylogenetic tree construction

1. Construct a phylogenetic tree from the filtered MSA.
2. Estimate branch support.
3. Identify major phylogenetic groups.
4. Examine the distribution of experimentally characterized proteins across the
tree.
5. Examine the distribution of Zn-binding site variants among phylogenetic
groups.

## Phylogenetic tree annotation

1. Annotate the phylogenetic tree with:
    - organism taxonomy;
    - COG assignment, where relevant;
    - presence of experimental structural data;
    - amino acid states at Zn-binding positions;
    - conservation or modification of the Zn-binding site.
2. Use taxonomic annotation to visualize the distribution of Zn-binding site
variants across major prokaryotic groups rather than treating organism taxonomy
itself as the reconstructed protein phylogeny.

## Evolutionary analysis of the Zn-binding site

1. Compare:
    - phylogenetic relationships;
    - amino acid states at Zn-binding positions;
    - Zn-binding site conservation patterns;
    - taxonomic distribution.
2. Determine:
    - the overall conservation of the Zn-binding site among Xni/FEN homologs;
    - which Zn-binding residues are highly conserved;
    - which positions exhibit substitutions;
    - which Zn-binding site variants occur;
    - whether particular variants are associated with specific phylogenetic
    groups.

## Final visualization and conclusions

1. Prepare:
    - the final filtered MSA with annotated functional positions;
    - an annotated phylogenetic tree;
    - a table of experimentally established Zn-binding residues;
    - representative structural views of the Zn-binding site;
    - a comparative visualization of Zn-binding site conservation across the
    phylogeny.
2. Summarize:
    - the degree of conservation of the Zn-binding site among Xni/FEN proteins
    within the analyzed COG0258 homologs;
    - the major amino acid substitutions affecting Zn-binding positions;
    - the phylogenetic distribution of Zn-binding site variants.
