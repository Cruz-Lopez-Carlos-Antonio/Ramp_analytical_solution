from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List

# ============================================================
# Configuration
# ============================================================
INPUT_FILE = Path("Corte_1.txt")
OUTPUT_FILE = Path("latex_tables.txt")

# The file can contain 6 tables; the paper needs the 3 reduced/even-time tables.
TABLES_TO_EXPORT = [4, 5, 6]

# APE formatting:
# The APE values are already formatted in the input file.
# Therefore, this script preserves the mantissa exactly as written:
#   7.20e-14   -> $7.20\times 10^{-14}$
#   1.7060e-16 -> $1.7060\times 10^{-16}$
#   1.80e-8    -> $1.80\times 10^{-8}$

# If True, the code only boldfaces the developed solution after the first
# differing character. If the numbers are identical as written, no bold is added.
BOLDFACE_FROM_FIRST_DIFFERENCE = True


@dataclass
class Row:
    t: str
    reference: str       # Runge--Kutta
    developed: str       # Proposed/developed solution
    ape: str


@dataclass
class TableBlock:
    number: int
    title: str
    rows: List[Row]


TEMPLATES: Dict[int, Dict[str, str]] = {
    4: {
        "placement": "h!",
        "caption": r"Neutron density values using \texttt{SciPy/NumPy} (16 digit precision), for the \textbf{Cold Start-up Scenario} and $a=0.0001\ \mathrm{s^{-1}}$.",
        "label": r"Table:n_density_SciPy",
        "quantity_header": r"Neutron density, $n(t)$, normalized",
        "developed_header": r"Developed solution (SciPy/NumPy)",
    },
    5: {
        "placement": "h!",
        "caption": r"Neutron density using \texttt{mpmath} (32 digit precision), for the \textbf{Cold Start-up Scenario} and $a=0.0001\ \mathrm{s^{-1}}$.",
        "label": r"Table:n_density_mpmath_reduced",
        "quantity_header": r"Neutron density, $n(t)$, normalized",
        "developed_header": r"Developed solution (\texttt{mpmath})",
    },
    6: {
        "placement": "t!",
        "caption": r"Delayed neutron's concentration using \texttt{SciPy/NumPy} (16 digit precision), for the \textbf{Cold Start-up Scenario} and $a=0.0001\ \mathrm{s^{-1}}$.",
        "label": r"Table:C_precursors_SciPy",
        "quantity_header": r"Precursor concentration, $C(t)$, normalized",
        "developed_header": r"Developed solution (SciPy/NumPy)",
    },
}


NUMBER_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
TABLE_TITLE_RE = re.compile(r"^Tabla\s+(\d+)\.\s*(.+?)\s*$")


def parse_tables(text: str) -> Dict[int, TableBlock]:
    """Parse all tables of the form used in comparative_tables.txt."""
    tables: Dict[int, TableBlock] = {}
    current_number: int | None = None
    current_title = ""
    current_rows: List[Row] = []

    def flush_current() -> None:
        nonlocal current_number, current_title, current_rows
        if current_number is not None:
            tables[current_number] = TableBlock(current_number, current_title, current_rows)
        current_number = None
        current_title = ""
        current_rows = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue

        title_match = TABLE_TITLE_RE.match(line)
        if title_match:
            flush_current()
            current_number = int(title_match.group(1))
            current_title = title_match.group(2)
            current_rows = []
            continue

        if current_number is None:
            continue

        parts = line.split()
        if len(parts) == 4 and all(NUMBER_RE.match(p) for p in parts):
            current_rows.append(Row(t=parts[0], reference=parts[1], developed=parts[2], ape=parts[3]))

    flush_current()
    return tables


def first_difference_index(a: str, b: str) -> int | None:
    """Return index where strings first differ; if one is a prefix, return prefix length."""
    limit = min(len(a), len(b))
    for i in range(limit):
        if a[i] != b[i]:
            return i
    if len(a) != len(b):
        return limit
    return None


def bold_from_first_difference(reference: str, developed: str) -> str:
    """Bold the part of developed beginning with its first difference from reference."""
    if not BOLDFACE_FROM_FIRST_DIFFERENCE:
        return developed

    idx = first_difference_index(reference, developed)
    if idx is None:
        return developed

    if idx >= len(developed):
        # The developed number is a strict prefix of the reference. There are
        # no visible remaining developed digits to boldface.
        return developed

    return developed[:idx] + r"\textbf{" + developed[idx:] + "}"


def latex_ape(value: str) -> str:
    """
    Convert decimal/scientific notation to a LaTeX math expression.

    Important:
    The mantissa is preserved exactly as it appears in the input file.
    This avoids losing significant trailing zeros, for example:
        1.7060e-16 -> $1.7060\times 10^{-16}$
        1.80e-8    -> $1.80\times 10^{-8}$
    """
    v = value.strip()

    try:
        if float(v) == 0.0:
            return r"$0$"
    except ValueError:
        return v

    if "e" in v.lower():
        mantissa_str, exponent_str = re.split(r"[eE]", v)
        exponent = int(exponent_str)
        return rf"${mantissa_str}\times 10^{{{exponent}}}$"

    return rf"${v}$"


def latex_row(row: Row) -> str:
    developed = bold_from_first_difference(row.reference, row.developed)
    ape = latex_ape(row.ape)
    return (
        f"{row.t:<2} & {row.reference} \n"
        f"   & {developed} \n"
        f"   & {ape} \\\\" 
    )


def render_table(table: TableBlock) -> str:
    cfg = TEMPLATES[table.number]
    rows = "\n".join(latex_row(row) + ("\n" if i != len(table.rows) - 1 else "") for i, row in enumerate(table.rows))
    return rf"""\begin{{table}}[{cfg['placement']}]
\centering
\small
\caption{{{cfg['caption']}}}
\label{{{cfg['label']}}}
\begin{{tabular}}{{c c c c}}
\hline
\makecell{{$t$ \\ {{[s]}}}} 
  & \multicolumn{{2}}{{c}}{{{cfg['quantity_header']}}} 
  & \makecell{{APE \\ {{[\,\%\,]}}}} \\
& Reference (Runge--Kutta 4th order) 
& {cfg['developed_header']} 
& \\
\hline
{rows}
\hline
\end{{tabular}}
\end{{table}}"""


def main() -> None:
    text = INPUT_FILE.read_text(encoding="utf-8")
    tables = parse_tables(text)

    missing = [n for n in TABLES_TO_EXPORT if n not in tables]
    if missing:
        raise RuntimeError(f"Missing expected table(s): {missing}")

    latex_tables = [render_table(tables[n]) for n in TABLES_TO_EXPORT]
    OUTPUT_FILE.write_text("\n\n".join(latex_tables) + "\n", encoding="utf-8")

    print(f"Created {OUTPUT_FILE} with tables {TABLES_TO_EXPORT}.")


if __name__ == "__main__":
    main()
