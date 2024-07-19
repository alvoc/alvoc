.PHONY: docs-start docs-build 

docs-start:
	@echo "Starting MkDocs server..."
	@export TERMYNAL_PREPROCESSOR_PRIORITY=26 && poetry run mkdocs serve

docs-build:
	@echo "Building MkDocs documentation..."
	@export TERMYNAL_PREPROCESSOR_PRIORITY=26 && poetry run mkdocs build
