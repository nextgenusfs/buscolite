# Contributing

Contributions to BUSCOlite are welcome! Here's how you can contribute:

## Setting Up Development Environment

1. Fork the repository on GitHub

2. Clone your fork locally:

   ```bash
   git clone https://github.com/your-username/buscolite.git
   cd buscolite
   ```

3. Install the package in development mode:

   ```bash
   pip install -e .
   ```

4. Install development dependencies:

   ```bash
   pip install pytest pytest-cov black flake8
   ```

## Running Tests

BUSCOlite uses pytest for testing. To run the tests:

```bash
# Run all tests
pytest

# Run with coverage report
pytest --cov=buscolite

# Run only unit tests
pytest tests/test_*.py

# Run only integration tests
pytest tests/integration/
```

## Code Style

BUSCOlite follows the PEP 8 style guide. You can use black and flake8 to ensure your code meets these standards:

```bash
# Format code with black
black buscolite tests

# Check code style with flake8
flake8 buscolite tests
```

## Pull Request Process

1. Create a new branch for your feature or bugfix:

   ```bash
   git checkout -b feature-name
   ```

2. Make your changes and commit them with descriptive commit messages

3. Push your branch to your fork:

   ```bash
   git push origin feature-name
   ```

4. Create a pull request on GitHub

5. Wait for the automated tests to pass

6. Address any review comments

## Adding Documentation

BUSCOlite uses MkDocs with Material theme for documentation. To build the documentation:

1. Install MkDocs and dependencies:

   ```bash
   pip install mkdocs-material mkdocstrings[python]
   ```

2. Serve the documentation locally:

   ```bash
   mkdocs serve
   ```

3. View the documentation in your browser at `http://127.0.0.1:8000`

4. Build the documentation:

   ```bash
   mkdocs build
   ```

## Adding New Features

When adding new features, please also add:

1. **Tests** for the new feature
2. **Documentation** for the new feature
3. **An entry** in the changelog

## Reporting Issues

If you find a bug or have a feature request, please open an issue on GitHub:

[https://github.com/nextgenusfs/buscolite/issues](https://github.com/nextgenusfs/buscolite/issues)

When reporting bugs, please include:

- Your operating system and version
- Python version
- BUSCOlite version
- Steps to reproduce the issue
- Expected behavior
- Actual behavior
- Any error messages or logs
