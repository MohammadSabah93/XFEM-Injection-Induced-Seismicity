# Contributing

Contributions that improve reproducibility, numerical robustness, documentation, or verification are welcome.

## Reporting an issue

Please include enough information to reproduce the problem:

- implementation branch (`Contact_Lagrange-Multiplier` or `Contact_Penalty`);
- MATLAB release and operating system;
- relevant model, friction, mesh, and time-stepping parameters;
- expected and observed behavior;
- complete error message or convergence information; and
- a minimal reproducible case when possible.

Please do not include confidential, proprietary, or personally identifiable data.

## Pull requests

Keep changes focused and explain their numerical or scientific motivation. For changes that affect model behavior, please include an appropriate verification, convergence comparison, or benchmark result.

A useful pull request should:

1. identify the affected formulation and files;
2. describe the mathematical or numerical change;
3. document any new parameters or assumptions;
4. preserve or explicitly update the documented default example; and
5. avoid committing generated result files, temporary MATLAB files, or large binary outputs.

Contributions are distributed under the repository's MIT License.
