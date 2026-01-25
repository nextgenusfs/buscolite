"""
SVG-based plotting for BUSCO assessment results.

This module provides simple, dependency-free plotting using only Python's
standard library xml.etree.ElementTree for SVG generation.
"""

import sys
import xml.etree.ElementTree as ET


class BuscoPlot:
    """Generate SVG bar chart for BUSCO assessment results."""

    # BUSCO color scheme
    COLORS = {
        "complete_single": "#56B4E9",  # Light blue
        "complete_duplicated": "#3492C7",  # Darker blue
        "fragmented": "#F0E442",  # Yellow
        "missing": "#F04442",  # Red
    }

    def __init__(self, width=800, height=400):
        """
        Initialize the BUSCO plot.

        Parameters
        ----------
        width : int
            Width of the SVG in pixels (default: 800)
        height : int
            Height of the SVG in pixels (default: 400)
        """
        self.width = width
        self.height = height
        self.margin = {"top": 60, "right": 40, "bottom": 80, "left": 120}
        self.plot_width = width - self.margin["left"] - self.margin["right"]
        self.plot_height = height - self.margin["top"] - self.margin["bottom"]

    def create_svg(self, stats, config, output_file):
        """
        Create an SVG bar chart from BUSCO statistics.

        Parameters
        ----------
        stats : dict
            Dictionary containing BUSCO statistics with keys:
            'total', 'single-copy', 'duplicated', 'fragmented', 'missing'
        config : dict
            Configuration dictionary with lineage information
        output_file : str
            Path to output SVG file
        """
        # Calculate percentages
        total = stats["total"]
        complete_single = stats["single-copy"]
        complete_dup = stats["duplicated"]
        fragmented = stats["fragmented"]
        missing = stats["missing"]

        complete_single_pct = (complete_single / total * 100) if total > 0 else 0
        complete_dup_pct = (complete_dup / total * 100) if total > 0 else 0
        fragmented_pct = (fragmented / total * 100) if total > 0 else 0
        missing_pct = (missing / total * 100) if total > 0 else 0

        # Create SVG root
        svg = ET.Element(
            "svg",
            {
                "width": str(self.width),
                "height": str(self.height),
                "xmlns": "http://www.w3.org/2000/svg",
                "version": "1.1",
            },
        )

        # Add white background
        ET.SubElement(
            svg,
            "rect",
            {
                "width": str(self.width),
                "height": str(self.height),
                "fill": "white",
            },
        )

        # Add title
        title_text = "BUSCO Assessment Results"
        title = ET.SubElement(
            svg,
            "text",
            {
                "x": str(self.width / 2),
                "y": "30",
                "text-anchor": "middle",
                "font-family": "Arial, sans-serif",
                "font-size": "20",
                "font-weight": "bold",
                "fill": "#333",
            },
        )
        title.text = title_text

        # Add subtitle with lineage info
        subtitle_text = f"Lineage: {config.get('name', 'Unknown')} (n={total})"
        subtitle = ET.SubElement(
            svg,
            "text",
            {
                "x": str(self.width / 2),
                "y": "50",
                "text-anchor": "middle",
                "font-family": "Arial, sans-serif",
                "font-size": "12",
                "fill": "#666",
            },
        )
        subtitle.text = subtitle_text

        # Calculate bar dimensions
        bar_height = 60
        bar_y = self.margin["top"] + (self.plot_height - bar_height) / 2

        # Draw stacked horizontal bar
        x_offset = self.margin["left"]

        # Complete single-copy
        if complete_single_pct > 0:
            width = (complete_single_pct / 100) * self.plot_width
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x_offset),
                    "y": str(bar_y),
                    "width": str(width),
                    "height": str(bar_height),
                    "fill": self.COLORS["complete_single"],
                },
            )
            x_offset += width

        # Complete duplicated
        if complete_dup_pct > 0:
            width = (complete_dup_pct / 100) * self.plot_width
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x_offset),
                    "y": str(bar_y),
                    "width": str(width),
                    "height": str(bar_height),
                    "fill": self.COLORS["complete_duplicated"],
                },
            )
            x_offset += width

        # Fragmented
        if fragmented_pct > 0:
            width = (fragmented_pct / 100) * self.plot_width
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x_offset),
                    "y": str(bar_y),
                    "width": str(width),
                    "height": str(bar_height),
                    "fill": self.COLORS["fragmented"],
                },
            )
            x_offset += width

        # Missing
        if missing_pct > 0:
            width = (missing_pct / 100) * self.plot_width
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x_offset),
                    "y": str(bar_y),
                    "width": str(width),
                    "height": str(bar_height),
                    "fill": self.COLORS["missing"],
                },
            )

        # Add percentage labels on the bar
        x_offset = self.margin["left"]
        label_y = bar_y + bar_height / 2 + 5

        # Helper function to add label if percentage is significant
        def add_label(pct, count, x_pos, width):
            if pct > 5:  # Only show label if segment is wide enough
                label_x = x_pos + width / 2
                label = ET.SubElement(
                    svg,
                    "text",
                    {
                        "x": str(label_x),
                        "y": str(label_y),
                        "text-anchor": "middle",
                        "font-family": "Arial, sans-serif",
                        "font-size": "12",
                        "font-weight": "bold",
                        "fill": "white",
                    },
                )
                label.text = f"{pct:.1f}%"

        if complete_single_pct > 0:
            width = (complete_single_pct / 100) * self.plot_width
            add_label(complete_single_pct, complete_single, x_offset, width)
            x_offset += width

        if complete_dup_pct > 0:
            width = (complete_dup_pct / 100) * self.plot_width
            add_label(complete_dup_pct, complete_dup, x_offset, width)
            x_offset += width

        if fragmented_pct > 0:
            width = (fragmented_pct / 100) * self.plot_width
            add_label(fragmented_pct, fragmented, x_offset, width)
            x_offset += width

        if missing_pct > 0:
            width = (missing_pct / 100) * self.plot_width
            add_label(missing_pct, missing, x_offset, width)

        # Add legend
        legend_y = bar_y + bar_height + 40
        legend_x = self.margin["left"]
        legend_spacing = 180

        categories = [
            ("Complete (S)", self.COLORS["complete_single"], complete_single, complete_single_pct),
            ("Complete (D)", self.COLORS["complete_duplicated"], complete_dup, complete_dup_pct),
            ("Fragmented (F)", self.COLORS["fragmented"], fragmented, fragmented_pct),
            ("Missing (M)", self.COLORS["missing"], missing, missing_pct),
        ]

        for i, (label, color, count, pct) in enumerate(categories):
            x = legend_x + (i * legend_spacing)

            # Color box
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x),
                    "y": str(legend_y),
                    "width": "15",
                    "height": "15",
                    "fill": color,
                },
            )

            # Label text
            text = ET.SubElement(
                svg,
                "text",
                {
                    "x": str(x + 20),
                    "y": str(legend_y + 12),
                    "font-family": "Arial, sans-serif",
                    "font-size": "11",
                    "fill": "#333",
                },
            )
            text.text = f"{label}: {count} ({pct:.1f}%)"

        # Add summary line
        summary_y = legend_y + 30
        complete_pct = complete_single_pct + complete_dup_pct
        summary_text = f"C:{complete_pct:.1f}%[S:{complete_single_pct:.1f}%,D:{complete_dup_pct:.1f}%],F:{fragmented_pct:.1f}%,M:{missing_pct:.1f}%,n:{total}"

        summary = ET.SubElement(
            svg,
            "text",
            {
                "x": str(self.width / 2),
                "y": str(summary_y),
                "text-anchor": "middle",
                "font-family": "monospace",
                "font-size": "11",
                "fill": "#333",
            },
        )
        summary.text = summary_text

        # Write to file
        tree = ET.ElementTree(svg)
        # ET.indent() was added in Python 3.9, so only use it if available
        if sys.version_info >= (3, 9):
            ET.indent(tree, space="  ")
        tree.write(output_file, encoding="utf-8", xml_declaration=True)


def generate_plot(results, missing, stats, config, output_file):
    """
    Generate a BUSCO assessment plot for a single result.

    Parameters
    ----------
    results : dict
        Dictionary of BUSCO results
    missing : list
        List of missing BUSCO IDs
    stats : dict
        Statistics dictionary
    config : dict
        Configuration dictionary
    output_file : str
        Path to output SVG file
    """
    plotter = BuscoPlot(width=800, height=400)
    plotter.create_svg(stats, config, output_file)


def generate_multi_plot(datasets, output_file):
    """
    Generate a comparative BUSCO assessment plot for multiple datasets.

    This function creates a stacked bar chart comparing BUSCO results across
    multiple samples, similar to BUSCO's official plotting functionality.

    Parameters
    ----------
    datasets : list of dict
        List of dictionaries, each containing:
        - 'name': str - Sample name
        - 'stats': dict - Statistics dictionary with keys:
          'total', 'single-copy', 'duplicated', 'fragmented', 'missing'
        - 'config': dict - Configuration dictionary (optional, uses first if not provided)
    output_file : str
        Path to output SVG file

    Example
    -------
    datasets = [
        {'name': 'Sample1', 'stats': {...}, 'config': {...}},
        {'name': 'Sample2', 'stats': {...}, 'config': {...}},
    ]
    generate_multi_plot(datasets, 'comparison.svg')
    """
    import xml.etree.ElementTree as ET

    if not datasets or len(datasets) == 0:
        raise ValueError("At least one dataset is required")

    # Calculate dimensions based on number of datasets
    n_datasets = len(datasets)
    bar_height = 40
    bar_spacing = 10
    margin = {"top": 80, "right": 40, "bottom": 100, "left": 200}
    width = 1000
    plot_height = n_datasets * (bar_height + bar_spacing)
    height = plot_height + margin["top"] + margin["bottom"]
    plot_width = width - margin["left"] - margin["right"]

    # BUSCO color scheme
    colors = {
        "complete_single": "#56B4E9",
        "complete_duplicated": "#3492C7",
        "fragmented": "#F0E442",
        "missing": "#F04442",
    }

    # Create SVG root
    svg = ET.Element(
        "svg",
        {
            "width": str(width),
            "height": str(height),
            "xmlns": "http://www.w3.org/2000/svg",
            "version": "1.1",
        },
    )

    # Add white background
    ET.SubElement(
        svg,
        "rect",
        {
            "width": str(width),
            "height": str(height),
            "fill": "white",
        },
    )

    # Add title
    title = ET.SubElement(
        svg,
        "text",
        {
            "x": str(width / 2),
            "y": "30",
            "text-anchor": "middle",
            "font-family": "Arial, sans-serif",
            "font-size": "20",
            "font-weight": "bold",
            "fill": "#333",
        },
    )
    title.text = "BUSCO Assessment Results - Comparison"

    # Add subtitle with lineage info (from first dataset)
    if "config" in datasets[0]:
        config = datasets[0]["config"]
        subtitle_text = f"Lineage: {config.get('name', 'Unknown')}"
        subtitle = ET.SubElement(
            svg,
            "text",
            {
                "x": str(width / 2),
                "y": "50",
                "text-anchor": "middle",
                "font-family": "Arial, sans-serif",
                "font-size": "12",
                "fill": "#666",
            },
        )
        subtitle.text = subtitle_text

    # Draw bars for each dataset
    for idx, dataset in enumerate(datasets):
        stats = dataset["stats"]
        name = dataset.get("name", f"Dataset {idx + 1}")

        total = stats["total"]
        complete_single = stats["single-copy"]
        complete_dup = stats["duplicated"]
        fragmented = stats["fragmented"]
        missing = stats["missing"]

        # Calculate percentages
        complete_single_pct = (complete_single / total * 100) if total > 0 else 0
        complete_dup_pct = (complete_dup / total * 100) if total > 0 else 0
        fragmented_pct = (fragmented / total * 100) if total > 0 else 0
        missing_pct = (missing / total * 100) if total > 0 else 0

        # Calculate bar position
        bar_y = margin["top"] + idx * (bar_height + bar_spacing)

        # Add dataset label
        label = ET.SubElement(
            svg,
            "text",
            {
                "x": str(margin["left"] - 10),
                "y": str(bar_y + bar_height / 2 + 5),
                "text-anchor": "end",
                "font-family": "Arial, sans-serif",
                "font-size": "12",
                "fill": "#333",
            },
        )
        label.text = name

        # Draw stacked bar
        x_offset = margin["left"]

        # Complete single-copy
        if complete_single_pct > 0:
            w = (complete_single_pct / 100) * plot_width
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x_offset),
                    "y": str(bar_y),
                    "width": str(w),
                    "height": str(bar_height),
                    "fill": colors["complete_single"],
                },
            )
            # Add percentage label if segment is wide enough
            if complete_single_pct > 5:
                text_elem = ET.SubElement(
                    svg,
                    "text",
                    {
                        "x": str(x_offset + w / 2),
                        "y": str(bar_y + bar_height / 2 + 4),
                        "text-anchor": "middle",
                        "font-family": "Arial, sans-serif",
                        "font-size": "10",
                        "font-weight": "bold",
                        "fill": "white",
                    },
                )
                text_elem.text = f"{complete_single_pct:.1f}%"
            x_offset += w

        # Complete duplicated
        if complete_dup_pct > 0:
            w = (complete_dup_pct / 100) * plot_width
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x_offset),
                    "y": str(bar_y),
                    "width": str(w),
                    "height": str(bar_height),
                    "fill": colors["complete_duplicated"],
                },
            )
            if complete_dup_pct > 5:
                text_elem = ET.SubElement(
                    svg,
                    "text",
                    {
                        "x": str(x_offset + w / 2),
                        "y": str(bar_y + bar_height / 2 + 4),
                        "text-anchor": "middle",
                        "font-family": "Arial, sans-serif",
                        "font-size": "10",
                        "font-weight": "bold",
                        "fill": "white",
                    },
                )
                text_elem.text = f"{complete_dup_pct:.1f}%"
            x_offset += w

        # Fragmented
        if fragmented_pct > 0:
            w = (fragmented_pct / 100) * plot_width
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x_offset),
                    "y": str(bar_y),
                    "width": str(w),
                    "height": str(bar_height),
                    "fill": colors["fragmented"],
                },
            )
            if fragmented_pct > 5:
                text_elem = ET.SubElement(
                    svg,
                    "text",
                    {
                        "x": str(x_offset + w / 2),
                        "y": str(bar_y + bar_height / 2 + 4),
                        "text-anchor": "middle",
                        "font-family": "Arial, sans-serif",
                        "font-size": "10",
                        "font-weight": "bold",
                        "fill": "white",
                    },
                )
                text_elem.text = f"{fragmented_pct:.1f}%"
            x_offset += w

        # Missing
        if missing_pct > 0:
            w = (missing_pct / 100) * plot_width
            ET.SubElement(
                svg,
                "rect",
                {
                    "x": str(x_offset),
                    "y": str(bar_y),
                    "width": str(w),
                    "height": str(bar_height),
                    "fill": colors["missing"],
                },
            )
            if missing_pct > 5:
                text_elem = ET.SubElement(
                    svg,
                    "text",
                    {
                        "x": str(x_offset + w / 2),
                        "y": str(bar_y + bar_height / 2 + 4),
                        "text-anchor": "middle",
                        "font-family": "Arial, sans-serif",
                        "font-size": "10",
                        "font-weight": "bold",
                        "fill": "white",
                    },
                )
                text_elem.text = f"{missing_pct:.1f}%"

    # Add legend at the bottom
    legend_y = margin["top"] + plot_height + 30
    legend_x = margin["left"]
    legend_spacing = 180

    categories = [
        ("Complete (S)", colors["complete_single"]),
        ("Complete (D)", colors["complete_duplicated"]),
        ("Fragmented (F)", colors["fragmented"]),
        ("Missing (M)", colors["missing"]),
    ]

    for i, (label, color) in enumerate(categories):
        x = legend_x + (i * legend_spacing)

        # Color box
        ET.SubElement(
            svg,
            "rect",
            {
                "x": str(x),
                "y": str(legend_y),
                "width": "15",
                "height": "15",
                "fill": color,
            },
        )

        # Label text
        text = ET.SubElement(
            svg,
            "text",
            {
                "x": str(x + 20),
                "y": str(legend_y + 12),
                "font-family": "Arial, sans-serif",
                "font-size": "11",
                "fill": "#333",
            },
        )
        text.text = label

    # Write to file
    tree = ET.ElementTree(svg)
    # ET.indent() was added in Python 3.9, so only use it if available
    if sys.version_info >= (3, 9):
        ET.indent(tree, space="  ")
    tree.write(output_file, encoding="utf-8", xml_declaration=True)
