#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess


def test_buscolite_help():
    """Test that the CLI help command works."""
    result = subprocess.run(
        ["buscolite", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "usage:" in result.stdout
    assert "buscolite" in result.stdout


def test_buscolite_version():
    """Test that the CLI version command works."""
    result = subprocess.run(
        ["buscolite", "--version"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "buscolite" in result.stdout
