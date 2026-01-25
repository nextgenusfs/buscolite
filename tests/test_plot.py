"""
Tests for the SVG plotting module.
"""

import os
import tempfile
import xml.etree.ElementTree as ET

import pytest

from buscolite.plot import BuscoPlot, generate_plot


def test_busco_plot_init():
    """Test BuscoPlot initialization."""
    plotter = BuscoPlot(width=800, height=400)
    assert plotter.width == 800
    assert plotter.height == 400
    assert plotter.plot_width == 800 - 120 - 40  # width - left - right margins
    assert plotter.plot_height == 400 - 60 - 80  # height - top - bottom margins


def test_create_svg():
    """Test SVG creation with sample data."""
    stats = {
        "total": 100,
        "single-copy": 80,
        "duplicated": 10,
        "fragmented": 5,
        "missing": 5,
    }
    config = {
        "name": "test_lineage_odb12",
        "creation_date": "2024-01-01",
        "number_of_species": "10",
        "number_of_BUSCOs": "100",
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
        output_file = f.name

    try:
        plotter = BuscoPlot(width=800, height=400)
        plotter.create_svg(stats, config, output_file)

        # Check that file was created
        assert os.path.exists(output_file)

        # Parse the SVG and check basic structure
        tree = ET.parse(output_file)
        root = tree.getroot()

        # Check it's an SVG
        assert root.tag.endswith("svg")
        assert root.attrib["width"] == "800"
        assert root.attrib["height"] == "400"

        # Check for rectangles (bars and legend boxes)
        rects = root.findall(".//{http://www.w3.org/2000/svg}rect")
        assert len(rects) > 0  # Should have background + bars + legend boxes

        # Check for text elements (title, labels, legend)
        texts = root.findall(".//{http://www.w3.org/2000/svg}text")
        assert len(texts) > 0

        # Check that title contains "BUSCO Assessment Results"
        title_found = False
        for text in texts:
            if text.text and "BUSCO Assessment Results" in text.text:
                title_found = True
                break
        assert title_found

    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)


def test_generate_plot():
    """Test the generate_plot helper function."""
    results = {
        "BUSCO1": {"status": "complete", "hit": "gene1", "bitscore": 100.0},
        "BUSCO2": {"status": "complete", "hit": "gene2", "bitscore": 90.0},
        "BUSCO3": {"status": "duplicated", "hit": "gene3", "bitscore": 85.0},
        "BUSCO4": {"status": "fragmented", "hit": "gene4", "bitscore": 70.0},
    }
    missing = ["BUSCO5"]
    stats = {
        "total": 5,
        "single-copy": 2,
        "duplicated": 1,
        "fragmented": 1,
        "missing": 1,
    }
    config = {
        "name": "test_lineage_odb12",
        "creation_date": "2024-01-01",
        "number_of_species": "10",
        "number_of_BUSCOs": "5",
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
        output_file = f.name

    try:
        generate_plot(results, missing, stats, config, output_file)

        # Check that file was created
        assert os.path.exists(output_file)

        # Parse and verify it's valid SVG
        tree = ET.parse(output_file)
        root = tree.getroot()
        assert root.tag.endswith("svg")

    finally:
        # Clean up
        if os.path.exists(output_file):
            os.remove(output_file)


def test_svg_with_zero_values():
    """Test SVG creation when some categories have zero values."""
    stats = {
        "total": 100,
        "single-copy": 100,
        "duplicated": 0,
        "fragmented": 0,
        "missing": 0,
    }
    config = {
        "name": "perfect_assembly",
        "creation_date": "2024-01-01",
        "number_of_species": "10",
        "number_of_BUSCOs": "100",
    }

    with tempfile.NamedTemporaryFile(mode="w", suffix=".svg", delete=False) as f:
        output_file = f.name

    try:
        plotter = BuscoPlot()
        plotter.create_svg(stats, config, output_file)

        # Should still create valid SVG
        assert os.path.exists(output_file)
        tree = ET.parse(output_file)
        root = tree.getroot()
        assert root.tag.endswith("svg")

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def test_svg_colors():
    """Test that the correct BUSCO colors are used."""
    assert BuscoPlot.COLORS["complete_single"] == "#56B4E9"
    assert BuscoPlot.COLORS["complete_duplicated"] == "#3492C7"
    assert BuscoPlot.COLORS["fragmented"] == "#F0E442"
    assert BuscoPlot.COLORS["missing"] == "#F04442"

