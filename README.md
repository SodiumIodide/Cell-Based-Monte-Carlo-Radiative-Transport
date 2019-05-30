# Cell Based Monte Carlo Radiative Transport

Deterministic and probability-based solutions to a time-dependent radiative transport model have been developed. Those solutions are attempted to be matched using a Monte Carlo particle-based method within a geometric cell. Future work hopefully includes coupling more than one cell via streaming.

## Three methods detailed

1. Emission is tracked as the particles travel within each time-step, these terms are coalesced into a total emission term

   **NOT WORKING**

2. Time-steps are split into discrete material, rather than tracked as transition of each particle

3. Emission is tracked at the end of the time-step in an explicit fashion
