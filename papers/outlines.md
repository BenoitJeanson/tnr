# Article 1
## Possible titles
- Risk-based approach to the Transmission Network Reconfiguration problem for congestion management.

# Outlines
## Introduction
- SOTA
    - The N-1 rule consists in mastering the consequences of contingencies rather than avoiding load loss at any cost.
    - Literature review on TNR applied to congestion management with MILP (cf. go-no Go).
    - + Worst case approach. No loss of connectivity is permitted.
- Contributions
    
## A risk-based approach to the TNR
- The root of the N-1 lays in the risk management.
- The common understanding of the rule is "no loss on a single contingency" is misleading:
    - The load is always connected to a final element that may trip
    - Ensuring total security, in all circumstances - especially during maintenance outages - would be very costly.
- The EU-SOGL does not require such a highly constrained standard
- Therefore, in some circumstances, the transmission grid is operated with pockets:
    - This happens mainly in "subtransmission" grids.
    - 3 levels operational connectedness: fully meshed (VHV) / with some pockets (HV) / radial (DSO)
- A pocket is rarely alone: the cascade of "pocketing" neighboring areas when weakening one link.
- Illustration of a cascade

## Modeling
Principles:
- DC-PF
- constraints: flow limit only
- cost function: $ risk = probability \times loss $

### (a PST-based model of the DC-PF)
if included, use Simon's paper.

### a branch-only MILP model
- DC-PF with branch-opening parameters
- N connectedness constraint
- N-1 connectedness management
    - identifying potential lost areas
    - balancing
- big-M implementation

### (a Bus-splitting MILP model)
- DC-PF with a library of possible predefined bus reconfiguration
- N connectedness constraint
- N-1 connectedness management
    - identifying potential lost areas
    - balancing
- big-M implementation

## Results
- Cases explained:
    - a N reconfigration required by a N-1 constraint (IEEE-14)
    - idem with necessity for pocketing (IEEE-14)
    - cascade of pockets (tbd)
- (OPF) vs Branch-only vs bus splitting
 
 ## Conclusion / discussion
 - WoW!! Awesome results!
 - Security limits breach is the red line the operators shall not cross. The security of supply comes next. Therefore, their main focus is the risk management: potentially putting at risk part of the supply to ensure being within the security limits. It shall be optimized first, before any other concern such as Joule losses.
 - Further developments:
    - scale
    - include other physical values such as voltage and short circuit.

# Pending questions
- use of PST model or phase (not yet implemented)
- OPF? 