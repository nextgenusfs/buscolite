"""
Tests for the log module.
"""

import os
import logging
import tempfile
from buscolite.log import startLogging


def test_startlogging_stdout():
    """Test startLogging function with stdout only."""
    # Call the function
    logger = startLogging()

    # Check that the logger is configured correctly
    assert isinstance(logger, logging.Logger)
    assert logger.level == logging.DEBUG

    # Check that there is one handler (stdout)
    assert len(logger.handlers) == 1
    assert isinstance(logger.handlers[0], logging.StreamHandler)
    assert logger.handlers[0].level == logging.INFO


def test_startlogging_with_file():
    """Test startLogging function with a log file."""
    # Create a temporary file for the log
    with tempfile.NamedTemporaryFile(suffix=".log", delete=False) as temp_log:
        log_path = temp_log.name

    try:
        # Call the function with a log file
        logger = startLogging(logfile=log_path)

        # Check that the logger is configured correctly
        assert isinstance(logger, logging.Logger)
        assert logger.level == logging.DEBUG

        # Check that there is at least one file handler
        # Note: The number of handlers may vary depending on previous test runs

        # Check the file handler
        file_handler = None
        for handler in logger.handlers:
            if isinstance(handler, logging.FileHandler):
                file_handler = handler
                break

        assert file_handler is not None
        assert file_handler.level == logging.DEBUG
        assert file_handler.baseFilename == log_path

        # Write a log message
        logger.info("Test log message")

        # Check that the message was written to the file
        with open(log_path, "r") as f:
            log_content = f.read()
            assert "Test log message" in log_content

    finally:
        # Clean up the temporary file
        if os.path.exists(log_path):
            os.remove(log_path)
