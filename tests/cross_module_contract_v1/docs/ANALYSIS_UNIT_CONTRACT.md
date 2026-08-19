# Analysis Unit Contract

Cells are the measurement layer. Animals/donors are the inference layer. Modules may enter one animal-level interpretation table, but they do not share one numeric scale.

```text
Kaile v1.1: single cells -> cell-level youth_score -> donor arithmetic mean for native deployment
Zihan v2.1: single-cell counts -> donor pseudobulk -> donor Youth Score
Kei Risk: single cells -> cell-level three-axis Risk scores -> animal-condition distribution summaries
Geneformer: single-cell embeddings -> cell-level diagnostic values -> donor-level validation required
Mashiro Risk: animal-level visual exists; exact pairing, threshold, CSV, and script remain TBD
```

`aggregate(score(cell)) != score(aggregate(expression))`. Kaile cell-score aggregation and Zihan pseudobulk scoring are different estimands.

