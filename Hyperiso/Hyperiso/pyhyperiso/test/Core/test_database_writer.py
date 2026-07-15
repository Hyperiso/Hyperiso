from pathlib import Path

import pytest

from pyhyperiso.core.Common.ParamId import ParamId
from pyhyperiso.core.Core.DatabaseWriter import DatabaseWriter


class _FakeCppWriter:
    def __init__(self):
        self.calls = []

    def write(self, destination):
        self.calls.append(("write", destination))

    def write_blocks(self, destination, block_names):
        self.calls.append(("write_blocks", destination, block_names))

    def write_parameters(self, destination, parameter_ids):
        self.calls.append(("write_parameters", destination, parameter_ids))


def _writer_with_fake_backend():
    writer = DatabaseWriter.__new__(DatabaseWriter)
    writer._cpp_obj = _FakeCppWriter()
    return writer


def test_database_writer_dispatches_supported_formats(tmp_path):
    writer = _writer_with_fake_backend()

    writer.write(tmp_path / "database.json")
    writer.write_blocks(tmp_path / "blocks.yaml", ["MASS", "SMINPUTS"])
    writer.write_parameters(
        tmp_path / "parameters.flha",
        [ParamId(block="MASS", code=25)],
    )

    assert writer._cpp_obj.calls[0] == (
        "write",
        str(tmp_path / "database.json"),
    )
    assert writer._cpp_obj.calls[1] == (
        "write_blocks",
        str(tmp_path / "blocks.yaml"),
        ["MASS", "SMINPUTS"],
    )
    assert writer._cpp_obj.calls[2][0:2] == (
        "write_parameters",
        str(tmp_path / "parameters.flha"),
    )
    assert len(writer._cpp_obj.calls[2][2]) == 1


@pytest.mark.parametrize("suffix", [".json", ".yaml", ".yml", ".lha", ".slha", ".flha"])
def test_database_writer_accepts_documented_suffixes(suffix):
    destination = DatabaseWriter._destination(Path("export").with_suffix(suffix))
    assert destination.endswith(suffix)


def test_database_writer_rejects_unsupported_suffix():
    writer = _writer_with_fake_backend()
    with pytest.raises(ValueError, match="Unsupported database export suffix"):
        writer.write("database.txt")


def test_database_writer_rejects_empty_selections():
    writer = _writer_with_fake_backend()
    with pytest.raises(ValueError, match="block_names"):
        writer.write_blocks("blocks.json", [])
    with pytest.raises(ValueError, match="parameter_ids"):
        writer.write_parameters("parameters.json", [])


def test_database_writer_rejects_non_paramid_entries():
    writer = _writer_with_fake_backend()
    with pytest.raises(TypeError, match="ParamId"):
        writer.write_parameters("parameters.json", [object()])


def test_database_writer_rejects_string_as_block_collection():
    writer = _writer_with_fake_backend()
    with pytest.raises(TypeError, match="iterable"):
        writer.write_blocks("blocks.json", "MASS")


def test_database_writer_rejects_non_string_block_names():
    writer = _writer_with_fake_backend()
    with pytest.raises(TypeError, match="only strings"):
        writer.write_blocks("blocks.json", ["MASS", object()])


def test_database_writer_native_json_integration(tmp_path, monkeypatch):
    """Initialize Core and verify that the native writer creates valid JSON."""
    import json

    from pyhyperiso.Common import Model
    from pyhyperiso.Core import HyperisoConfig, HyperisoMaster

    monkeypatch.setenv("HYPERISO_CACHE_ROOT", str(tmp_path / "cache"))
    master = HyperisoMaster()
    master.init("lha/testInput.flha", HyperisoConfig(model=Model.SM))

    destination = tmp_path / "database.json"
    DatabaseWriter().write(destination)

    payload = json.loads(destination.read_text(encoding="utf-8"))
    assert payload
    assert any(str(key).upper() == "MASS" for key in payload)
