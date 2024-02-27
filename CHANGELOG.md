### Version 0.1.5 (February 27, 2024) ###
Added a new feature that evaluates how representative a set of taxa is! The representatives for evaluation need to be indicated with the `--prior` regex and they can be evaluated using the `--evaluate` output option.

### Version 0.1.4 (August 10, 2023) ###
Added an option to save the partitioning of a tree by PARNAS to a file (see `--clusters`).

### Version 0.1.3 (July 9, 2023) ###
Added support for Python 3.11 and other improvements.

### Version 0.1.2 (December 22, 2022) ###
- Bug fix: Closed [Issue #1](https://github.com/flu-crew/parnas/issues/1).
- Added safe guards for taxon weights: negative weights or weights exceeding 1000 are no longer allowed.
