# Testing BUSCOlite

This document describes the testing setup for BUSCOlite and how to run the tests.

## Testing Framework

The BUSCOlite project uses pytest for testing. The tests are organized as follows:

- `tests/test_*.py`: Unit tests for individual functions and modules
- `tests/integration/`: Integration tests for the entire workflow
  - `tests/integration/test_buscolite_cli.py`: Tests for the command-line interface
  - `tests/integration/test_buscolite_workflow.py`: Tests for the workflow using the Python API
  - `tests/integration/test_buscolite_api.py`: Tests for the Python API functions
  - `tests/integration/test_buscolite_end_to_end.py`: End-to-end tests using subprocess

## Running Tests

To run the tests, you need to install pytest and the package in development mode:

```bash
# Set up the development environment (installs pytest, pre-commit, and other tools)
./scripts/setup_dev.sh

# Or manually install the development dependencies
pip install -e ".[dev]"

# Run all tests
pytest

# Run with coverage report
pytest --cov=buscolite

# Generate HTML coverage report
python scripts/run_coverage.py

# Run only unit tests
pytest tests/test_*.py

# Run only integration tests
pytest tests/integration/
```

## Test Data

The test data is stored in the following directories:

- `tests/integration/data/`: Test data for integration tests
  - `test_genome.fasta`: A small test genome
  - `test_proteins.fasta`: A small test proteome
  - `mock_lineage/`: A mock BUSCO lineage directory

## Requirements

The tests require the following dependencies:

- pytest
- pytest-cov
- buscolite (installed in development mode)
- augustus (for genome mode tests)
- miniprot (for genome mode tests)

Some tests are skipped if the required dependencies are not installed.

## Current Test Coverage

The current test coverage is as follows:

- `__init__.py`: 100%
- `log.py`: 100%
- `fasta.py`: 91%
- `fastx.py`: 81%
- `gff.py`: 73%
- `help_formatter.py`: 64%
- `__main__.py`: 63%
- `utilities.py`: 30%
- `busco.py`: 18%
- `search.py`: 16%
- `augustus.py`: 11%

Overall coverage: 46%

## Continuous Integration

The tests are run automatically on GitHub Actions when changes are pushed to the repository. The workflow is defined in the `.github/workflows/tests.yml` file.

The GitHub Actions workflow does the following:

1. Runs on push and pull requests to the main/master branch
2. Tests with multiple Python versions (3.8, 3.9, 3.10, 3.11)
3. Installs the required dependencies (including miniprot)
4. Runs the unit tests and generates a coverage report
5. Runs the integration tests that don't require external dependencies
6. Uploads the coverage report to Codecov
7. Runs linting checks with flake8 and black

You can see the status of the tests on the [GitHub Actions page](https://github.com/nextgenusfs/buscolite/actions) and the coverage report on the [Codecov page](https://codecov.io/gh/nextgenusfs/buscolite).

## Code Formatting and Linting

The project uses pre-commit hooks to enforce code formatting and linting standards. The following tools are used:

- **black**: For code formatting with a line length of 100 characters
- **isort**: For sorting imports
- **flake8**: For linting

To run the formatting and linting checks manually:

```bash
# Run all pre-commit hooks
pre-commit run --all-files

# Run specific hooks
pre-commit run black --all-files
pre-commit run isort --all-files
pre-commit run flake8 --all-files
```

## Adding New Tests

When adding new features or fixing bugs, please also add tests to ensure the code works correctly. Follow these guidelines:

1. For unit tests, add a test function to the appropriate test file in the `tests/` directory.
2. For integration tests, add a test function to the appropriate test file in the `tests/integration/` directory.
3. Make sure the test function has a descriptive name and docstring.
4. Use pytest fixtures to set up test data and clean up after tests.
5. Use pytest.mark.skipif to skip tests that require external dependencies that might not be installed.
6. Make sure your code passes all pre-commit checks before submitting a pull request.

## Running Tests with Specific Python Versions

To run the tests with a specific Python version, you can use a virtual environment:

```bash
# Create a virtual environment with Python 3.9
python3.9 -m venv venv-py39
source venv-py39/bin/activate

# Install dependencies
pip install pytest pytest-cov
pip install -e .

# Run tests
pytest
```
