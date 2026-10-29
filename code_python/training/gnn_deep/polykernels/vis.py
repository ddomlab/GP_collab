"""Plotting and visualization helpers for the optimization pipeline."""

import io
from pathlib import Path
import re

from loguru import logger
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import rdkit.Chem
import rdkit.Chem.Draw
import rdkit.Chem.rdDepictor
import seaborn as sns  # type: ignore

# Output resolution for every raster we write. Anything rasterized on its way into a
# figure (RDKit molecule drawings) must be generated at this density, not upscaled to it.
FIG_DPI = 300

# Canvas width RDKit's default font and bond-width settings are tuned for.
_MOL_RASTER_BASE_WIDTH = 500

# Vertical space reserved for a horizontal colourbar strip, in inches.
_CBAR_HEIGHT_IN = 0.55

# RDKit highlights every bond whose two endpoints are both highlighted -- which, since
# `draw_molecule_with_scores` highlights every atom, means every bond -- drawing it at
# `bondLineWidth * highlightBondWidthMultiplier`. At RDKit's default multiplier of 8 that
# band is wide enough to merge adjacent atoms' discs into a continuous ribbon, so the
# per-atom colour that carries the attribution is what gets lost. Kept just above the
# bond line itself, so a bond reads as an outlined stroke rather than a coloured tube.
_HIGHLIGHT_BOND_WIDTH_MULTIPLIER = 2

# Highlight disc radius in RDKit's bond-length units: a floor every atom gets, plus a
# term scaled by |score| / max|score|. Both sit above RDKit's 0.3 default so the disc,
# not the bond band, is the dominant mark. The ceiling stays under 0.5 -- half a bond
# length -- so neighbouring discs stop short of swallowing the bond between them.
_HIGHLIGHT_RADIUS_BASE = 0.32
_HIGHLIGHT_RADIUS_SPAN = 0.16

# Atom-label and atom-index type size. RDKit caps label height at `maxFontSize` = 40 px by
# default, which is tuned for a ~500 px canvas; on the 2400 px raster these figures use,
# that cap binds long before the label is proportionate to the molecule, so heteroatom
# symbols and indices come out tiny. Lifting the cap hands sizing back to `baseFontSize`,
# which is a fraction of *bond length* and so tracks the drawing at any canvas size --
# verified to stay legible from a 21-atom monomer up to dataset1's 186-atom pair.
_ATOM_FONT_SIZE_CAP = -1  # -1 disables the cap
_ATOM_FONT_SIZE = 0.60
# Index annotations, as a fraction of the atom-label size. Above ~0.75 RDKit starts
# dropping indices it can no longer place, which loses information rather than adding it.
_ATOM_INDEX_FONT_SCALE = 0.70

# Ink and neutral tones shared by the figure helpers below.
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#b5b3aa"

# Hard ceiling on a figure dimension. A heatmap with hundreds of rows silently grows to
# tens of inches, which any viewer then scales down until the tick labels are unreadable;
# callers should subset their rows rather than rely on this.
_MAX_FIG_IN = 24.0


def set_visualization_style() -> None:
    """Set global plotting style parameters for polished scientific figures."""
    sns.set_theme(
        style="ticks",
        palette="husl",
        rc={
            "axes.grid": True,
            "grid.color": "#e5e5e5",
            "grid.linestyle": "--",
            "grid.linewidth": 0.6,
            "savefig.dpi": FIG_DPI,
            "savefig.transparent": True,
            "savefig.pad_inches": 0.1,
        },
    )


def save_figure(
    fig: plt.Figure, path: Path, *, vector_overlay: str | None = None
) -> None:
    """Save a matplotlib figure as both PNG (solid white bg) and SVG (transparent).

    Args:
        vector_overlay: SVG markup for the figure's single raster image, if it has one.
            Matplotlib can only embed a bitmap into an axes, so a figure built around an
            `imshow` writes that bitmap into the SVG too and the "vector" output is
            vector in name only. When given, the bitmap is swapped for this markup in the
            SVG so it stays resolution-independent; the PNG is unaffected.
    """
    # PNG with opaque white background
    fig.savefig(
        path,
        dpi=FIG_DPI,
        bbox_inches="tight",
        facecolor="white",
        transparent=False,
    )
    # SVG with transparent background
    svg_path = path.with_suffix(".svg")
    if vector_overlay is None:
        fig.savefig(svg_path, bbox_inches="tight", transparent=True)
    else:
        buf = io.BytesIO()
        fig.savefig(buf, format="svg", bbox_inches="tight", transparent=True)
        svg_path.write_text(
            _replace_svg_raster(buf.getvalue().decode("utf-8"), vector_overlay),
            encoding="utf-8",
        )
    plt.close(fig)
    logger.info(f"Saved figure to {path} and {svg_path}")


# Artists tagged with this gid are the ones `vector_overlay` replaces. Matplotlib writes
# a gid straight through to the SVG element's id, which is the only reliable way to pick
# the intended raster out of a figure -- a colourbar is also drawn as an <image>.
_VECTOR_SLOT_GID = "vector-overlay-slot"

# Matplotlib emits an imshow raster as a single self-closing <image> whose placement is
# expressed in flipped coordinates, e.g.
#   <image xlink:href="data:..." transform="scale(1 -1) translate(0 -156.24)"
#          x="7.2" y="-7.2" width="223.2" height="156.24"/>
_MPL_FLIP_TRANSFORM = re.compile(
    r"transform=\"scale\(1 -1\) translate\(0 (-?[0-9.]+)\)\""
)


def _replace_svg_raster(svg: str, overlay: str) -> str:
    """Swap the `_VECTOR_SLOT_GID` raster <image> in a matplotlib SVG for vector markup.

    A nested `<svg>` carrying the overlay's own viewBox re-establishes its coordinate
    system, so the drawing scales into the slot the bitmap occupied without any need to
    reconcile the two coordinate spaces by hand.
    """
    tags = [
        t for t in re.findall(r"<image\b[^>]*/>", svg)
        if f'id="{_VECTOR_SLOT_GID}"' in t
    ]
    if len(tags) != 1:
        raise ValueError(
            f"Expected exactly one <image> tagged {_VECTOR_SLOT_GID}, found {len(tags)}. "
            f"vector_overlay requires the artist it replaces to carry that gid."
        )
    tag = tags[0]
    # Drop the data URI before reading attributes: base64 padding can end in `x=`, which
    # an attribute regex would otherwise mistake for the image's own x attribute.
    attrs = re.sub(r"xlink:href=\"[^\"]*\"", "", tag)

    flip = _MPL_FLIP_TRANSFORM.search(attrs)
    if not flip:
        raise ValueError(f"Unrecognized transform on matplotlib <image>: {attrs}")

    def attr(name: str) -> float:
        m = re.search(rf'\b{name}="(-?[0-9.]+)"', attrs)
        if not m:
            raise ValueError(f"No {name} attribute on matplotlib <image>: {attrs}")
        return float(m.group(1))

    x, y = attr("x"), attr("y")
    w, h = attr("width"), attr("height")
    # Under scale(1 -1) translate(0 -T) a local point (x, y) lands at (x, T - y), so the
    # top edge of the rect in root coordinates is T - y - height.
    top = -float(flip.group(1)) - y - h

    root = re.search(r"<svg\b[^>]*>", overlay)
    if not root:
        raise ValueError("vector_overlay contains no <svg> element.")
    open_tag = root.group(0)
    view_box = re.search(r"viewBox=['\"]([^'\"]*)['\"]", open_tag)
    if not view_box:
        raise ValueError("vector_overlay's <svg> has no viewBox, so it cannot be scaled.")

    # Rebuild the open tag with the slot's geometry; the overlay's own width/height would
    # otherwise fight the placement we just computed.
    stripped = re.sub(r"\b(width|height|x|y)=['\"][^'\"]*['\"]", "", open_tag)
    placed = (
        f'{stripped[:-1].rstrip()} x="{x}" y="{top}" width="{w}" height="{h}" '
        f'preserveAspectRatio="xMidYMid meet">'
    )
    return svg.replace(tag, placed + overlay[root.end():], 1)


def draw_molecule_with_scores(
    smiles: str,
    scores: np.ndarray,
    *,
    cmap_name: str,
    signed: bool,
    title: str,
    width: int = 2400,
    height: int = 1680,
) -> tuple[plt.Figure, str]:
    """Draw a molecule with each atom coloured by its attribution score.

    Ported from `scripts/Che1148_CodeLab8_IntegratedGradients.ipynb`, returning the figure
    instead of calling `plt.show()` so callers can hand it to `save_figure`.

    The molecule is drawn twice from identical options: a Cairo raster that matplotlib can
    place in an axes, and SVG markup that only `save_figure` can use. Pass the second
    return value as that function's `vector_overlay` to get a genuinely scalable SVG.

    Args:
        signed: Colour on a diverging scale centred at zero and size each highlight by
            |score|. When False, scores are treated as magnitudes on a one-sided scale.
        width: Pixel width of the RDKit raster. This is also the pixel width of the saved
            PNG, so `imshow` maps the bitmap roughly 1:1 instead of upscaling it. RDKit
            has no vector path into a matplotlib axes, so the only way to keep the
            molecule sharp at `FIG_DPI` is to rasterize it at the output size.

    Returns:
        The figure, and the molecule panel as standalone SVG markup.
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    rdkit.Chem.rdDepictor.Compute2DCoords(mol)

    scores = np.asarray(scores, dtype=float)
    if mol.GetNumAtoms() != len(scores):
        raise ValueError(
            f"{mol.GetNumAtoms()} atoms but {len(scores)} scores for SMILES: {smiles}"
        )

    if signed:
        vmax = float(np.max(np.abs(scores)) + 1e-12)
        norm = mpl.colors.TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
    else:
        vmax = float(np.max(scores) + 1e-12)
        norm = mpl.colors.Normalize(vmin=0.0, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)

    atom_colors, atom_radii = {}, {}
    for i, score in enumerate(scores):
        atom_colors[i] = tuple(float(v) for v in cmap(float(norm(float(score)))))
        atom_radii[i] = float(
            _HIGHLIGHT_RADIUS_BASE
            + _HIGHLIGHT_RADIUS_SPAN * abs(float(score)) / vmax
        )

    def render(drawer) -> None:
        opts = drawer.drawOptions()
        opts.addAtomIndices = True
        # RDKit's bond width is an absolute pixel count tuned for a canvas around
        # `_MOL_RASTER_BASE_WIDTH`; on a print-resolution canvas bonds come out hairline.
        opts.bondLineWidth = max(1, round(2 * width / _MOL_RASTER_BASE_WIDTH))
        opts.highlightBondWidthMultiplier = _HIGHLIGHT_BOND_WIDTH_MULTIPLIER
        opts.maxFontSize = _ATOM_FONT_SIZE_CAP
        opts.baseFontSize = _ATOM_FONT_SIZE
        opts.annotationFontScale = _ATOM_INDEX_FONT_SCALE
        # The SVG is composited into a transparent figure, so RDKit must not paint its
        # own opaque backdrop over the panel. Harmless for the raster, which matplotlib
        # draws over an explicitly white background anyway.
        opts.clearBackground = False
        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(range(mol.GetNumAtoms())),
            highlightAtomColors=atom_colors,
            highlightAtomRadii=atom_radii,
        )
        drawer.FinishDrawing()

    raster = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DCairo(width, height)
    render(raster)
    img = Image.open(io.BytesIO(raster.GetDrawingText()))

    # Same geometry and options, so the two outputs are the same drawing.
    vector = rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    render(vector)
    mol_svg = vector.GetDrawingText()

    # Size the axes from the raster so the molecule is not resampled on the way out.
    img_w_in, img_h_in = width / FIG_DPI, height / FIG_DPI
    fig, (ax_img, ax_cbar) = plt.subplots(
        2, 1,
        figsize=(img_w_in, img_h_in + _CBAR_HEIGHT_IN),
        gridspec_kw={"height_ratios": [img_h_in, _CBAR_HEIGHT_IN]},
        dpi=FIG_DPI,
    )
    ax_img.imshow(img, interpolation="antialiased").set_gid(_VECTOR_SLOT_GID)
    ax_img.axis("off")
    ax_img.grid(visible=False)
    ax_img.set_title(title, fontsize=10)

    mappable = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    mappable.set_array([])
    cbar = fig.colorbar(mappable, cax=ax_cbar, orientation="horizontal")
    cbar.set_label(
        "Signed node attribution" if signed else "Magnitude of node attribution"
    )
    fig.tight_layout()
    return fig, mol_svg


def plot_interval_bars(
    values: np.ndarray,
    errors: np.ndarray,
    labels: list[str],
    *,
    title: str,
    xlabel: str,
    secondary: np.ndarray | None = None,
    secondary_label: str = "",
) -> plt.Figure:
    """Horizontal bars with symmetric error bars, sorted by magnitude.

    Built for stochastic Shapley values, where the interval is the point of the plot rather
    than decoration: a feature whose credible interval straddles zero is one the model
    cannot place, which a bare bar chart would hide.

    Args:
        errors: half-width of the interval per bar (e.g. 1.96 * sd for 95%).
        secondary: optional second series drawn as open markers on the same rows, for
            contrasting two summaries of the same features.
    """
    set_visualization_style()
    order = np.argsort(np.abs(values))
    values, errors = np.asarray(values)[order], np.asarray(errors)[order]
    labels = [labels[i] for i in order]

    height = max(2.5, 0.32 * len(labels) + 1.4)
    fig, ax = plt.subplots(figsize=(8, height))
    y = np.arange(len(labels))
    # Diverging fill so sign is readable without consulting the axis.
    colors = ["#2a78d6" if v >= 0 else "#eb6834" for v in values]
    ax.barh(y, values, height=0.62, color=colors, zorder=3)
    ax.errorbar(
        values, y, xerr=errors, fmt="none", ecolor=INK_2, elinewidth=1.1, capsize=3,
        zorder=4,
    )
    if secondary is not None:
        ax.scatter(
            np.asarray(secondary)[order], y, facecolors="none", edgecolors=INK,
            s=32, linewidth=1.1, zorder=5, label=secondary_label,
        )
        ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax.axvline(0, color=MUTED, lw=0.9, zorder=2)
    ax.set_yticks(y, labels, fontsize=8)
    ax.set_xlabel(xlabel, fontsize=8.5)
    ax.set_title(title, loc="left", fontsize=11, color=INK, pad=10)
    ax.grid(axis="x", color=MUTED, alpha=0.4, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()
    return fig


def plot_attribution_heatmap(
    matrix: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    *,
    title: str,
    diverging: bool = True,
    cbar_label: str = "Attribution",
    annot: bool | np.ndarray = False,
    annot_fmt: str = ".2f",
) -> plt.Figure:
    """Heatmap of an attribution matrix (rows x columns).

    Args:
        diverging: Centre the colour scale at zero, for signed attributions. Set False for
            magnitudes, which are one-sided and would waste half a diverging ramp.
        annot: Print each cell's value inside it. Left off by default because the per-atom
            and per-feature heatmaps run to hundreds of cells, where the numbers overlap
            into noise; switch it on for small summary grids that are read as a table.
            An array of strings the shape of `matrix` prints those instead, e.g. "0.41 ± 0.06".
    """
    set_visualization_style()
    matrix = np.asarray(matrix, dtype=float)
    # A row label may carry a second line of context (dataset size, and so on). Scale the
    # band height with the tallest label so those lines get their own room instead of
    # overprinting the neighbouring row.
    label_lines = max((str(lab).count("\n") + 1 for lab in row_labels), default=1)
    height = min(_MAX_FIG_IN, max(2.5, 0.22 * label_lines * len(row_labels) + 1.5))
    width = min(_MAX_FIG_IN, max(6.0, 0.32 * len(col_labels) + 3.0))
    text_annot = isinstance(annot, np.ndarray)
    if text_annot:
        # A "0.412 ± 0.031" cell needs several times a bare value's width, and the row
        # labels beside it take their own share, so budget both: the longest annotation at
        # 8pt per column, the longest row-label line at 7pt, and the colour bar.
        cell_in = 0.075 * max(len(str(a)) for a in annot.ravel())
        label_in = 0.055 * max(len(line) for lab in row_labels for line in str(lab).split("\n"))
        width = min(_MAX_FIG_IN, max(width, cell_in * len(col_labels) + label_in + 1.5))
    fig, ax = plt.subplots(figsize=(width, height))

    if diverging:
        vmax = float(np.nanmax(np.abs(matrix))) or 1.0
        kwargs = {"cmap": "RdBu_r", "vmin": -vmax, "vmax": vmax, "center": 0.0}
    else:
        kwargs = {"cmap": "viridis", "vmin": 0.0}

    sns.heatmap(
        matrix,
        ax=ax,
        xticklabels=col_labels,
        yticklabels=row_labels,
        cbar_kws={"label": cbar_label},
        linewidths=0.2,
        linecolor="white",
        annot=annot,
        fmt="" if text_annot else annot_fmt,
        # Scale with the cells so the labels stay inside them on a wide grid.
        annot_kws={"fontsize": 8} if text_annot or annot else None,
        **kwargs,
    )
    ax.set_title(title, fontsize=11, loc="left")
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    ax.tick_params(axis="y", rotation=0, labelsize=7)
    ax.grid(visible=False)
    fig.tight_layout()
    return fig
