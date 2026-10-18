"""Execute the results notebook and export an HTML copy with embedded outputs."""
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/vegandiet-matplotlib")
os.environ.setdefault("JUPYTER_RUNTIME_DIR", "/tmp/vegandiet-jupyter-runtime")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks/figure1b_results.ipynb"


def main():
    notebook = nbformat.read(NOTEBOOK, as_version=4)
    nbformat.validate(notebook)
    client = NotebookClient(
        notebook,
        timeout=300,
        kernel_name=notebook.metadata.kernelspec.name,
        resources={"metadata": {"path": str(ROOT)}},
        allow_errors=False,
    )
    client.execute()
    nbformat.validate(notebook)
    if any(output.output_type == "error"
           for cell in notebook.cells if cell.cell_type == "code"
           for output in cell.outputs):
        raise RuntimeError("Notebook contains an execution error")
    temporary = NOTEBOOK.with_suffix(".ipynb.tmp")
    nbformat.write(notebook, temporary)
    temporary.replace(NOTEBOOK)
    exporter = HTMLExporter()
    exporter.exclude_input_prompt = True
    exporter.exclude_output_prompt = True
    html, _ = exporter.from_notebook_node(notebook)
    html_path = NOTEBOOK.with_suffix(".html")
    temporary = html_path.with_suffix(".html.tmp")
    temporary.write_text(html)
    temporary.replace(html_path)
    executed = sum(cell.cell_type == "code" for cell in notebook.cells)
    print(f"Executed {executed} code cells without errors.")
    print(f"Notebook: {NOTEBOOK}")
    print(f"HTML: {html_path}")


if __name__ == "__main__":
    main()
