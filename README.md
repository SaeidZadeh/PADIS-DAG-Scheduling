# PADIS: Precedence-Aware Deadline-Impact Scheduling

R implementation of the **PADIS** algorithm from:

> Alirezazadeh, S., Alibabaei, K., Rodrigues, N., & Pereira, A.
> *Precedence-Aware Deadline-Impact Scheduling (PADIS) for Time-Critical
> DAG Tasks on Multi-Core Multiprocessors.*
> Preprint, 2025.

## Overview

PADIS is a lightweight online heuristic for scheduling time-critical
DAG tasks on multiple multi-core processors. It maximizes the number of
accepted tasks (tasks that complete before their individual deadlines)
using a lexicographic priority key that combines:

- Task deadline (primary)
- Successor-release impact term *O_v* (tie-break)

Two variants are provided:
- **PADIS-ZC** — zero inter-processor communication (used for baseline comparison)
- **PADIS-CA** — communication-aware (see paper Section 4)

## Requirements

- R ≥ 4.0
- CRAN packages: `plot3D`, `vioplot`, `matrixStats`
- Bioconductor packages: `pcalg`, `Rgraphviz`

Install all dependencies by running:

```r
source("padis.R")   # auto-installs on first run
```

## License

MIT — see [LICENSE](LICENSE).

## Contact

Saeid Alirezazadeh — saeid.alirezazadeh@gmail.com
