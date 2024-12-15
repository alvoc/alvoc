import logging
from pathlib import Path
import typer
from functools import wraps
from rich.console import Console

# Common options
virus: str = typer.Argument(
    ..., help="Taxonomic ID of the virus, or path to the GenBank file"
)
outdir: Path = typer.Option(
    ".", help="Output directory for results and intermediate data"
)
logger = logging.getLogger("alvoc")

# Initialize Rich console
console = Console()

# Spinner decorator with dynamic text updates
# Spinner decorator with dynamic text updates and exclusive spinner output
def with_spinner(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        with console.status("[bold blue]Initializing...[/bold blue]", spinner="dots", spinner_style="blue") as status:
            # Custom logging handler to update spinner dynamically
            class SpinnerHandler(logging.Handler):
                def emit(self, record):
                    log_entry = self.format(record)
                    status.update(f"[bold blue]{log_entry}[/bold blue]")

            # Save existing handlers and temporarily disable them
            root_logger = logging.getLogger()  # Root logger
            existing_handlers = root_logger.handlers[:]
            root_logger.handlers = []  # Clear all existing handlers

            # Add spinner handler
            spinner_handler = SpinnerHandler()
            spinner_handler.setFormatter(logging.Formatter("%(message)s"))
            root_logger.addHandler(spinner_handler)

            try:
                result = func(*args, **kwargs)
            finally:
                # Restore original handlers
                root_logger.handlers = existing_handlers

            return result

    return wrapper

