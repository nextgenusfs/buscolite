import os
import shutil
import subprocess
import sys


def _get_buscolite_cmd():
    if shutil.which("buscolite"):
        return ["buscolite"]
    # Fallback to standard local script locations or module execution
    user_bin = os.path.expanduser("~/Library/Python/3.9/bin/buscolite")
    if os.path.exists(user_bin):
        return [user_bin]
    return [sys.executable, "-m", "buscolite"]


def test_buscolite_help():
    """Test that the CLI help command works."""
    result = subprocess.run(
        _get_buscolite_cmd() + ["--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "usage:" in result.stdout
    assert "buscolite" in result.stdout


def test_buscolite_version():
    """Test that the CLI version command works."""
    result = subprocess.run(
        _get_buscolite_cmd() + ["--version"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "buscolite" in result.stdout
