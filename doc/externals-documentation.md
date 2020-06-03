This document tracks external software fixes done in order to run it with the BLR pipeline.

# NAIBR

## Solved issues 

Date: 2020-06-03
- `mpmath` and `future` is not in original git environment description but needed to run.
- Needs to be run with work directory = NAIBR folder (Otherwise won't find its own modules).
