import logging


def init_logger(level=logging.DEBUG):
    logging.basicConfig(
        level=level,
        format="%(levelname)s - %(asctime)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[],
    )


def get_logger(name=None):
    return logging.getLogger(name) if name else logging.getLogger()
