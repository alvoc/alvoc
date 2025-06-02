# Developer Guide

This guide covers everything you need to know about developing alvoc.

## Project Structure

The project is organized as follows:

```
alvoc/
├── alvoc/              # Main package source code
│   ├── cli/           # Command line interface
│   ├── core/          # Core functionality
│   │   ├── amplicons/    # Amplicon analysis
│   │   ├── constellations/  # Constellation management
│   │   ├── utils/          # Shared utilities
│   │   └── variants/       # Variant analysis
│   └── main.py         # Entry point
├── docs/              # Documentation
├── examples/          # Example scripts
└── tests/             # Test suite
```

### Key Components

- **CLI Module** (`alvoc/cli/`): Handles command-line interface using [Typer](https://typer.tiangolo.com/). 
- **Core Module** (`alvoc/core/`): 
  - `amplicons/`: Amplicon metrics and visualization
  - `constellations/`: Mutation constellation management
  - `variants/`: Variant/lineage analysis
  - `utils/`: Shared utilities like logging and file handling

## Development Setup

1. Clone the repository:
```bash
git clone https://github.com/alvoc/alvoc
cd alvoc
```

2. Install development dependencies:
```bash
uv pip install -e ".[dev]"
```

## Running Tests

We use pytest for testing. The test suite is in the `tests/` directory.

To run tests:

```bash
pytest
```

Key testing guidelines:
- Place tests in the `tests/` directory
- Name test files with `test_` prefix
- Use fixtures for reusable test components
- Mock external dependencies using `unittest.mock`

Example of a good test:

```python
@pytest.fixture
def mock_seq_record():
    features = [
        SeqFeature(
            FeatureLocation(start=0, end=10),
            type="gene",
            qualifiers={"gene": ["gene1"]},
        )
    ]
    return SeqRecord(Seq("ATGC"), id="test_id", features=features)

def test_extract_gene_info(mock_seq_record):
    expected = {"gene1": [0, 10]}
    result = extract_gene_info(mock_seq_record)
    assert result == expected
```

## Documentation

Documentation is built using MkDocs with the Material theme. If you have UV installed, you can use our makefile to streamline this process

### Local Development

1. Start the documentation server:
```bash
make docs-start
```

2. View at http://127.0.0.1:8000

### Writing Documentation

- Documentation source files are in `docs/`
- Use Markdown format
- API reference is auto-generated from docstrings
- Add new pages to `mkdocs.yml` under `nav:`
 
## Release Process

Whenever you push a change to main, our CI/CD pipelines will automatically parse through commits and bump the version as required. Note that currently, alvoc will not increment past the major version number 0. To disable that behaviour, you can edit the pyproject.toml file:

```toml
major_version_zero = true # -> set this to false
```

The release goes through the followings steps:
 
1. **Testing**: 
   - Runs on pull requests and pushes to main
   - Executes test suite
   - Checks code formatting

2. **Documentation**:
   - Builds and deploys documentation on merge to main
   - Deploys to GitHub Pages

3. **Release**:
   - Triggered by new version tags
   - Builds and publishes package to PyPI
   - Builds and publishes Docker image to GHCR
   - Updates documentation

## Docker Support

The project includes multi-stage Dockerfile:

1. Base stage with Python environment
2. Build stage for dependencies
3. Runtime stage for minimal production image

Build the Docker image:
```bash
docker build -t alvoc .
```

Run the container:
```bash
docker run -it --rm alvoc
```

## Best Practices

1. **Code Style**:
   - Follow PEP 8
   - Use type hints
   - Document functions with docstrings

2. **Git Workflow**:
   - Use feature branches
   - Write meaningful commit messages following Conventional Commits
   - Keep commits atomic and focused

3. **Testing**:
   - Write tests for new features
   - Maintain good test coverage
   - Mock external dependencies

4. **Documentation**:
   - Keep API docs up-to-date
   - Document complex procedures
   - Include examples