import logging


def init_logger(level=logging.DEBUG, log_file="app.log"):
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler()],
    )


def get_logger(name=None):
    return logging.getLogger(name) if name else logging.getLogger()
