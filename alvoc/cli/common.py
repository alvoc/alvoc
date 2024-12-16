import logging
from pathlib import Path
import typer
from functools import wraps
from rich.console import Console
import time

# Common options
virus: str = typer.Argument(
    ..., help="Path to a GenBank file or NCBI Entrez query terms"
)
outdir: Path = typer.Option(
    ".", "--outdir", "-o", help="Output directory for results and intermediate data"
)

# Initialize Rich console
console = Console()


# Spinner decorator with dynamic text updates and exclusive spinner output
def with_spinner(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        with console.status(
            "[bold blue]Initializing...[/bold blue]",
            spinner="dots",
            spinner_style="blue",
        ) as status:
            # Custom logging handler to update spinner dynamically
            class SpinnerHandler(logging.Handler):
                def emit(self, record):
                    log_entry = self.format(record)
                    status.update(f"[bold blue]{log_entry}[/bold blue]")

            # Add spinner handler
            root_logger = logging.getLogger()
            spinner_handler = SpinnerHandler()
            spinner_handler.setFormatter(logging.Formatter("%(message)s"))
            spinner_handler.setLevel(logging.INFO)  # Only log INFO and above
            root_logger.addHandler(spinner_handler)

            try:
                result = func(*args, **kwargs)
            finally:
                time.sleep(0.5)  # Ensure spinner is visible
                status.update("[bold green]Task completed![/bold green]")
                root_logger.removeHandler(spinner_handler)  # Clean up handler

            return result

    return wrapper
