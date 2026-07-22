# GTEP Progressive Hedging Example

This directory contains an example configuration and local development driver for
running Progressive Hedging (PH) on the GTEP model using representative periods
as scenarios.

## Overview

The PH formulation treats each representative period as one scenario.

For example, with four representative periods:

```text
Scenario 1 = representative period 1
Scenario 2 = representative period 2
Scenario 3 = representative period 3
Scenario 4 = representative period 4

python ra_to_gtep.py --egret-json /Users/jkskolf/Downloads/2030_pcm_case/WI_2030_sens_medloadg_baseexp_schedret_tradload_tx_base.json --meta-json /Users/jkskolf/Downloads/2030_pcm_case/WI_2030_ts_medloadg_baseexp_schedret_tradload_tx_base.json --h5 /Users/jkskolf/Downloads/2030_pcm_case/WI_2030_ts_medloadg_baseexp_schedret_tradload_tx_base.h5 --outdir initial_case_2030_test



python -m gtep.algorithms.progressive_hedging.torc_orchestrator \
  --config examples/progressive_hedging/gtep_ph_config.yaml \
  --lineage gtep_ph_123_bus_resil_week

  SQLX_OFFLINE=true cargo build --release --features "server-bin,mcp-server,dash,slurm-runner"

