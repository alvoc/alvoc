.PHONY: docs-start docs-build docs-deploy

docs-start:
	@echo "Starting MkDocs server..."
	@export TERMYNAL_PREPROCESSOR_PRIORITY=26 && uv run mkdocs serve

docs-build:
	@echo "Building MkDocs documentation..."
	@export TERMYNAL_PREPROCESSOR_PRIORITY=26 && uv run mkdocs build

docs-deploy:
	@echo "Building MkDocs documentation..."
	@export TERMYNAL_PREPROCESSOR_PRIORITY=26 && uv run mkdocs gh-deploy --force
