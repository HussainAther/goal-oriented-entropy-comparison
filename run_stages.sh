#!/bin/sh
# Independent stage execution, suitable for CI jobs with per-command time limits.
set -eu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
case "${1:-}" in
 perturbation) (cd code/perturbation && python run_ent03p.py);;
 features) cp code/perturbation/cases.csv code/features/inputs/cases.csv; (cd code/features && python run_ent03q.py);;
 confirmation) cp code/features/joined_cases.csv code/confirmation/discovery_cases.csv; (cd code/confirmation && python run_ent03v.py);;
 comparator) cp code/features/joined_cases.csv code/comparator/dependencies/discovery_cases.csv; cp code/confirmation/confirmation_cases.csv code/comparator/dependencies/confirmation_cases.csv; (cd code/comparator && python run_ent04b.py);;
 figures) python build_release.py;;
 *) echo 'Run stages in this order: perturbation, features, confirmation, comparator, figures' >&2; exit 2;;
esac
