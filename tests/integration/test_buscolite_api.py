#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import buscolite


def test_buscolite_import():
    """Test that buscolite can be imported."""
    assert buscolite.__version__ is not None


def test_buscolite_logger():
    """Test that the logger can be initialized."""
    from buscolite.log import startLogging

    logger = startLogging()
    assert logger is not None
