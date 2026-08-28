# RGA prediction report -- r570_subsample

*Generated 2026-08-28 22:00:30 UTC by `rgas_prediction.py` v0.0.1 (config v0.0.1).*

## 1. What this report shows

4 proteins were examined and 1 (25.00%) carry at least one feature associated with plant immune receptors. These are *candidates* identified from protein domains and topology: they are not experimentally validated resistance genes.

## 2. How the call was made

Six independent annotation tools are harmonised into a single controlled vocabulary of protein features (NB-ARC, TIR, RPW8, CC, LRR, kinase, LysM, transmembrane helix, signal peptide). Protein domains come from InterProScan and are matched by accession only, never by description text, because descriptions change between releases and match unrelated entries. Overlapping hits reported by several signature databases for the same region are merged before anything is counted, so one LRR seen by Pfam, SMART and Gene3D counts once. Transmembrane helices are taken from Phobius and DeepTMHMM, signal peptides from SignalP 6.0 and Phobius, and coiled coils from three channels: a domain-level profile HMM plus the DeepCoil2 and InterProScan Coils predictors. A helix predicted inside the signal peptide is discarded, because signal peptides are routinely mistaken for transmembrane helices. Each protein is then passed through an ordered list of mutually exclusive rules and receives the first class that fits, together with a written justification citing the exact signatures behind the call. Subcellular localisation from DeepLoc 2.0 never decides a class; it only raises or lowers the reported confidence and flags inconsistencies.

## 3. Run metadata

### Reproduce this run

The exact command that produced this report, quoted as it was invoked, so it can be pasted back into a shell from the repository root:

```bash
uv run python code/rgas/rgas_prediction.py \
    --interproscan r570_subsample_interpro.tsv \
    --phobius r570_subsample_phobius.tsv \
    --deeptmhmm r570_subsample_deeptmhmm.gff3 \
    --signalp r570_subsample_signalp6_predictions.txt \
    --deeploc r570_subsample_deeploc2.csv \
    --deepcoil deepcoil_data \
    --outdir rga_out/ \
    --organism-name r570_subsample \
    --workers 6
```

### Settings

- Output directory: `rga_out`
- Consensus policies: TM `union`, SP `signalp`, CC `union`
- Coiled-coil calling: threshold 0.5, minimum length 21 residues, maximum gap 2 residues
- Minimum LRR copies: 1

### Input files

| tool | path | available | size_bytes | n_lines | sha256 |
| --- | --- | --- | --- | --- | --- |
| interproscan | r570_subsample_interpro.tsv | True | 5008 | 31 | 95bed53cd0c71e4ac4c1e318a2b0046a43ebaec58c9e6cb2e0e9ffb7fdabc344 |
| phobius | r570_subsample_phobius.tsv | True | 237 | 5 | ae3d00876262b6a2a6130086743bcd8e32d9f8e1a945a770511c5eac334be127 |
| deeptmhmm | r570_subsample_deeptmhmm.gff3 | True | 1081 | 25 | 635653586030b35f8d02047502e86bd536dce9d9381a83e9d7ebcc39525d461a |
| signalp | r570_subsample_signalp6_predictions.txt | True | 1299 | 6 | 49c48c18f06698e2615eb301c4be379564d64bdfdcb755ac068fcb94588722ca |
| deeploc | r570_subsample_deeploc2.csv | True | 1679 | 5 | 2534ee272eaeed0b86eced99feb43a8a0e6536e91cb98da2479eecc64bd859d3 |
| deepcoil | deepcoil_data | True | 41630 | 4 | ca36d80d27cf1c454890f76e2196c36f54ab848a8ba1554162502834b934350b |

### Evidence channels

| channel | available |
| --- | --- |
| interproscan | True |
| phobius | True |
| deeptmhmm | True |
| signalp | True |
| deeploc | True |
| deepcoil | True |

## 4. Rules applied

| priority | rule_id | family | subclass | requires | requires_one_of | forbids | description |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | CNL | NLR | CNL | NB-ARC;CC;LRR | - | TIR;RPW8 | Coiled-coil NLR: CC + NB-ARC + LRR |
| 2 | TNL | NLR | TNL | NB-ARC;TIR;LRR | - | RPW8 | TIR NLR: TIR + NB-ARC + LRR |
| 3 | RNL | NLR | RNL | NB-ARC;RPW8;LRR | - | - | Helper NLR (ADR1/NRG1 type): RPW8 + NB-ARC + LRR |
| 4 | NL | NLR | NL | NB-ARC;LRR | - | CC;TIR;RPW8 | NB-ARC + LRR without an N-terminal CC/TIR/RPW8 domain |
| 5 | CN | NLR | CN | NB-ARC;CC | - | LRR;TIR;RPW8 | CC + NB-ARC, LRR not detected |
| 6 | TN | NLR | TN | NB-ARC;TIR | - | LRR;RPW8 | TIR + NB-ARC, LRR not detected |
| 7 | RN | NLR | RN | NB-ARC;RPW8 | - | LRR | RPW8 + NB-ARC, LRR not detected |
| 8 | N | NLR | N | NB-ARC | - | CC;TIR;RPW8;LRR | NB-ARC only |
| 9 | TX | NLR-associated | TX | TIR | - | NB-ARC;LRR | TIR-X / TIR-only: TIR without NB-ARC and without LRR |
| 10 | RX | NLR-associated | RX | RPW8 | - | NB-ARC;TIR | RPW8-X: RPW8 without NB-ARC |
| 11 | LRR-RLK | RLK | LRR-RLK | STTK;LRR | (TM OR SP) | NB-ARC;TIR;RPW8 | Kinase + LRR ectodomain + (TM or SP), no NB-ARC |
| 12 | LysM-RLK | RLK | LysM-RLK | STTK;LysM | (TM OR SP) | NB-ARC;TIR;RPW8;LRR | Kinase + LysM ectodomain + (TM or SP), no NB-ARC |
| 14 | LRR-RLP | RLP | LRR-RLP | LRR | (TM OR SP) | NB-ARC;TIR;RPW8;STTK | LRR ectodomain + (TM or SP), no kinase, no NB-ARC |
| 15 | LysM-RLP | RLP | LysM-RLP | LysM | (TM OR SP) | NB-ARC;TIR;RPW8;STTK;LRR | LysM ectodomain + (TM or SP), no kinase, no NB-ARC |
| 16 | other-RLP | RLP | other-RLP | - | (TM OR SP) AND (LRR OR LysM) | NB-ARC;TIR;RPW8;STTK;LRR;LysM | Non-LRR/non-LysM ectodomain + (TM or SP), no kinase, no NB-ARC |
| 17 | TM-CC | TM-CC | TM-CC | TM;CC | - | NB-ARC;STTK;LRR;LysM;TIR;RPW8 | Transmembrane + coiled coil, no NB-ARC/kinase/LRR/LysM |
| 18 | Other | Other | Other | - | - | - | Carries at least one core immune feature but fits no rule above |
| 19 | Non-RGA | Non-RGA | NA | - | - | - | No core immune feature detected |

## 5. Counts

### By family

| rga_family | n_proteins | percent_of_proteome |
| --- | --- | --- |
| Non-RGA | 3 | 75.0 |
| NLR | 1 | 25.0 |

### By subclass

| rga_family | rga_subclass | n_proteins | percent_of_proteome |
| --- | --- | --- | --- |
| NLR | CNL | 1 | 25.0 |
| Non-RGA | NA | 3 | 75.0 |

### Confidence of RGA calls

| confidence | n_proteins |
| --- | --- |
| medium | 1 |


### Confidence by subclass

| rga_subclass | high | medium | low | n_proteins |
| --- | --- | --- | --- | --- |
| CNL | 0 | 1 | 0 | 1 |

### Most frequent domain architectures among RGAs

| domain_architecture | n_proteins |
| --- | --- |
| CC-NB-ARC-LRR | 1 |

## 6. Coiled-coil evidence

The coiled coil is the least reliable feature in every published RGA pipeline, and it is the one that decides CNL against NL. Three channels are used here, and they are not of equal weight. The leading one is a curated profile HMM for a named domain (the Rx N-terminal domain, PF18052 / IPR041118), which carries the same kind of evidence as the NB-ARC model every NLR call already rests on. The other two, DeepCoil2 and InterProScan Coils, are biophysical propensity predictors: neither publishes a recommended score cut-off, and Simm et al. (2021), benchmarking coiled-coil predictors against the whole PDB, found a 30-fold spread in how many coiled coils they call and agreement with structure close to random. They are kept because they cover proteins no domain model reaches, and a call resting on them alone is graded down rather than hidden. The tables below show how much the channels disagree and how far the subclass counts move with the policy.

### DeepCoil2 versus InterProScan Coils (whole proteome)

| InterProScan Coils | DeepCoil2 | Rx domain | n_proteins |
| --- | --- | --- | --- |
| CC called | CC called | -- | 1 |
| no CC | CC called | -- | 0 |
| CC called | no CC | -- | 1 |
| no CC | no CC | -- | 2 |
| -- | -- | CC called | 1 |
| no CC | no CC | CC called | 0 |

### Subclass counts under each `--cc-policy`

| rga_subclass | rx_domain | deepcoil | coils | union | intersection |
| --- | --- | --- | --- | --- | --- |
| CNL | 1 | 0 | 1 | 1 | 0 |
| NA | 3 | 3 | 3 | 3 | 3 |
| NL | 0 | 1 | 0 | 0 | 1 |

### Sensitivity to the segment-calling parameters

| threshold | min_length | n_proteins_with_cc | n_segments |
| --- | --- | --- | --- |
| 0.2 | 14 | 2 | 4 |
| 0.2 | 21 | 1 | 1 |
| 0.2 | 28 | 1 | 1 |
| 0.5 | 14 | 1 | 1 |
| 0.5 | 21 | 1 | 1 |
| 0.5 | 28 | 1 | 1 |

## 7. Identifier reconciliation

| tool | n_ids | n_shared_with_proteome | n_absent_from_tool | n_not_in_proteome |
| --- | --- | --- | --- | --- |
| deepcoil | 4 | 4 | 0 | 0 |
| deeploc | 4 | 4 | 0 | 0 |
| deeptmhmm | 4 | 4 | 0 | 0 |
| interproscan | 4 | 4 | 0 | 0 |
| phobius | 4 | 4 | 0 | 0 |
| signalp | 4 | 4 | 0 | 0 |

## 8. Warnings

| warning | n_proteins |
| --- | --- |
| DeepLoc localisation (Cell membrane) is inconsistent with class CNL | 1 |

## 9. Top 1 RGA candidates

| protein_id | rga_family | rga_subclass | domain_architecture | n_lrr | predicted_localization | confidence |
| --- | --- | --- | --- | --- | --- | --- |
| SoffiXsponR570.10Eg022100.1.p | NLR | CNL | CC-NB-ARC-LRR | 1 | Cell membrane | medium |

## 10. References

1. Rody HVS, Bombardelli RGH, Creste S, Camargo LEA, Van Sluys M-A, Monteiro-Vitorello CB (2019). Genome survey of resistance gene analogs in sugarcane: genomic features and differential expression of the innate immune system from a smut-resistant genotype. BMC Genomics 20:809. doi:10.1186/s12864-019-6207-y
2. Li P, Quan X, Jia G, Xiao J, Cloutier S, You FM (2016). RGAugury: a pipeline for genome-wide prediction of resistance gene analogs (RGAs) in plants. BMC Genomics 17:852. doi:10.1186/s12864-016-3197-x
3. Sekhwal MK, Li P, Lam I, Wang X, Cloutier S, You FM (2015). Disease resistance gene analogs (RGAs) in plants. Int J Mol Sci 16:19248-19290. doi:10.3390/ijms160819248
4. Kourelis J, Sakai T, Adachi H, Kamoun S (2021). RefPlantNLR is a comprehensive collection of experimentally validated plant disease resistance proteins from the NLR family. PLoS Biology 19(10):e3001124. doi:10.1371/journal.pbio.3001124
5. Smith M, Jones JT, Hein I (2025). Resistify: a novel NLR classifier that reveals Helitron-associated NLR expansion in Solanaceae. Bioinform Biol Insights 19:11779322241308944. doi:10.1177/11779322241308944
6. Shiu S-H, Bleecker AB (2003). Expansion of the receptor-like kinase/Pelle gene family and receptor-like proteins in Arabidopsis. Plant Physiol 132:530-543. doi:10.1104/pp.103.021964
7. Jones JDG, Dangl JL (2006). The plant immune system. Nature 444:323-329. doi:10.1038/nature05286
8. Jones P et al. (2014). InterProScan 5: genome-scale protein function classification. Bioinformatics 30:1236-1240. doi:10.1093/bioinformatics/btu031
9. Blum M et al. (2025). InterPro: the protein sequence classification resource in 2025. Nucleic Acids Res 53:D444-D456. doi:10.1093/nar/gkae1082
10. Paysan-Lafosse T et al. (2025). The Pfam protein families database: embracing AI/ML. Nucleic Acids Res 53:D523-D534. doi:10.1093/nar/gkae997
11. Kall L, Krogh A, Sonnhammer ELL (2004). A combined transmembrane topology and signal peptide prediction method. J Mol Biol 338:1027-1036. doi:10.1016/j.jmb.2004.03.016
12. Hallgren J et al. (2022). DeepTMHMM predicts alpha and beta transmembrane proteins using deep neural networks. bioRxiv. doi:10.1101/2022.04.08.487609
13. Teufel F et al. (2022). SignalP 6.0 predicts all five types of signal peptides using protein language models. Nat Biotechnol 40:1023-1025. doi:10.1038/s41587-021-01156-3
14. Thumuluri V et al. (2022). DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. Nucleic Acids Res 50:W228-W234. doi:10.1093/nar/gkac278
15. Ludwiczak J, Winski A, Szczepaniak K, Alva V, Dunin-Horkawicz S (2019). DeepCoil - a fast and accurate prediction of coiled-coil domains in protein sequences. Bioinformatics 35(16):2790-2795. doi:10.1093/bioinformatics/bty1062
16. Lupas A, Van Dyke M, Stock J (1991). Predicting coiled coils from protein sequences. Science 252:1162-1164. doi:10.1126/science.252.5009.1162
17. Simm D, Hatje K, Waack S, Kollmar M (2021). Critical assessment of coiled-coil predictions based on protein structure data. Scientific Reports 11:12439. doi:10.1038/s41598-021-91886-w
