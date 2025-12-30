"""Common utility functions for pangenomium commands"""
import os
import sys
from loguru import logger

def create_directory(path):
    """Create directory if it doesn't exist"""
    os.makedirs(path, exist_ok=True)
    return path

def setup_directories(output_directory, subdirectory=None):
    """Setup standard directory structure for pangenomium workflows
    
    Args:
        output_directory: Base output directory
        subdirectory: Optional subdirectory name for command-specific isolation
                     (e.g., 'genome_clustering', 'protein_clustering')
    
    Returns dict with keys: project, output, intermediate, log, tmp
    
    Note: If subdirectory is provided, intermediate and tmp are namespaced
          but output and log remain shared for easier access to results
    """
    directories = {
        "project": create_directory(output_directory),
    }
    
    # Shared directories (all commands write here)
    directories["output"] = create_directory(os.path.join(directories["project"], "output"))
    directories["log"] = create_directory(os.path.join(directories["project"], "log"))
    
    # Command-specific directories (isolated per command)
    if subdirectory:
        directories["intermediate"] = create_directory(
            os.path.join(directories["project"], "intermediate", subdirectory)
        )
        directories["tmp"] = create_directory(
            os.path.join(directories["project"], "tmp", subdirectory)
        )
    else:
        directories["intermediate"] = create_directory(
            os.path.join(directories["project"], "intermediate")
        )
        directories["tmp"] = create_directory(
            os.path.join(directories["project"], "tmp")
        )
    
    # Set TMPDIR environment variable
    os.environ["TMPDIR"] = directories["tmp"]
    
    return directories

def setup_logger(log_directory=None, log_filename="pangenomium.log"):
    """Setup loguru logger with file and console outputs
    
    Args:
        log_directory: Directory for log file (if None, logs to console only)
        log_filename: Name of log file
    """
    # Remove default handler
    logger.remove()
    
    # Add console handler with INFO level
    logger.add(
        sys.stderr,
        format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
        level="INFO"
    )
    
    # Add file handler if log directory specified
    if log_directory:
        log_path = os.path.join(log_directory, log_filename)
        logger.add(
            log_path,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
            level="DEBUG"
        )
    
    return logger

def print_header(version, n_jobs, additional_info=None):
    """Print standard configuration header using loguru
    
    Args:
        version: Script version
        n_jobs: Number of threads
        additional_info: Dict of additional key-value pairs to print
    """
    from pyexeggutor import get_timestamp
    
    logger.info("="*80)
    logger.info("Configuration")
    logger.info("-"*80)
    logger.info(f"Python version: {sys.version.split()[0]}")
    logger.info(f"Python path: {sys.executable}")
    logger.info(f"Script version: {version}")
    logger.info(f"Moment: {get_timestamp()}")
    logger.info(f"Directory: {os.getcwd()}")
    logger.info(f"Commands: {' '.join(sys.argv)}")
    logger.info(f"Threads: {n_jobs}")
    
    if additional_info:
        for key, value in additional_info.items():
            logger.info(f"{key}: {value}")
    
    logger.info("="*80)
