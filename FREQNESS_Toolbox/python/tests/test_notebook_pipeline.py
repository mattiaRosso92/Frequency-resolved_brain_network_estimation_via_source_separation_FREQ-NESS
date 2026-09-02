import json
from pathlib import Path

from freqness import ALL_ANALYSES


def test_notebook_is_clean_and_covers_every_pipeline_analysis() -> None:
    notebook_path = (
        Path(__file__).resolve().parents[1] / "FREQNESS_NotebookPipeline.ipynb"
    )
    notebook = json.loads(notebook_path.read_text(encoding="utf-8"))
    cells = notebook["cells"]
    source = "\n".join("".join(cell.get("source", [])) for cell in cells)

    for analysis in ALL_ANALYSES:
        assert analysis in source

    assert "FREQNESS_InducedResponses" in source
    assert "single-condition companion" in source
    assert "unattended processing of every condition" in source

    for cell in cells:
        if cell["cell_type"] != "code":
            continue
        assert cell["execution_count"] is None
        assert cell["outputs"] == []
        compile("".join(cell["source"]), str(notebook_path), "exec")


def test_batch_and_package_pipeline_files_are_retained() -> None:
    python_root = Path(__file__).resolve().parents[1]
    batch_pipeline = python_root / "FREQNESS_MainPipeline.py"
    package_pipeline = (
        python_root / "src" / "freqness" / "FREQNESS_MainPipeline.py"
    )

    assert batch_pipeline.is_file()
    assert package_pipeline.is_file()
    assert "run_configured_pipeline" in batch_pipeline.read_text(encoding="utf-8")
    assert "def FREQNESS_MainPipeline" in package_pipeline.read_text(
        encoding="utf-8"
    )
