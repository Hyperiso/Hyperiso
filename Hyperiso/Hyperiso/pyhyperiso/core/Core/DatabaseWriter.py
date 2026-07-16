"""Export the initialized HyperIso Core database to structured files."""

from __future__ import annotations

import os
from collections.abc import Iterable
from pathlib import Path
from typing import Union

from pyhyperiso.phyperiso.pyhyperiso.core import DatabaseWriter as _CppDatabaseWriter
from pyhyperiso.core.Common.ParamId import ParamId

PathLike = Union[str, os.PathLike[str]]


class DatabaseWriter:
    """Write current Core parameters as JSON, YAML or Les Houches files.

    The writer operates on the database initialized by :class:`HyperisoMaster`.
    The destination suffix selects the serializer:

    - ``.json`` for JSON;
    - ``.yaml`` or ``.yml`` for YAML;
    - ``.lha``, ``.slha`` or ``.flha`` for Les Houches formats.

    JSON and YAML preserve an optional ``imaginary_value`` field for complex
    parameters. LHA-family exports preserve complex entries through ``IM...``
    companion blocks, unless an explicit companion block already exists in the
    database.

    Examples:
        >>> writer = DatabaseWriter()
        >>> writer.write("database.json")
        >>> writer.write_blocks("masses.yaml", ["MASS", "SMINPUTS"])
        >>> writer.write_parameters("inputs.slha", [ParamId(block="MASS", code=25)])
    """

    _SUPPORTED_SUFFIXES = {".json", ".yaml", ".yml", ".lha", ".slha", ".flha"}

    def __init__(self) -> None:
        """Create a writer bound to the current C++ Core database."""
        self._cpp_obj = _CppDatabaseWriter()

    @classmethod
    def _destination(cls, destination: PathLike) -> str:
        """Validate and normalize an output path for the C++ binding."""
        path = Path(os.fspath(destination)).expanduser()
        suffix = path.suffix.lower()
        if suffix not in cls._SUPPORTED_SUFFIXES:
            supported = ", ".join(sorted(cls._SUPPORTED_SUFFIXES))
            raise ValueError(
                f"Unsupported database export suffix {suffix!r}; expected one of {supported}."
            )
        return str(path)

    def write(self, destination: PathLike) -> None:
        """Export the complete initialized database.

        Args:
            destination: Output filename. Its suffix selects JSON, YAML or LHA.

        Raises:
            RuntimeError: If HyperIso has not been initialized or writing fails.
            ValueError: If the filename suffix is unsupported.
        """
        self._cpp_obj.write(self._destination(destination))

    def write_blocks(self, destination: PathLike, block_names: Iterable[str]) -> None:
        """Export only selected blocks.

        Args:
            destination: Output filename.
            block_names: Non-empty iterable of block names or aliases.

        Raises:
            TypeError: If block names are not supplied as strings.
            ValueError: If no block is supplied, the suffix is unsupported, or a
                requested block does not exist.
            RuntimeError: If HyperIso is not initialized or writing fails.
        """
        if isinstance(block_names, (str, bytes)):
            raise TypeError("block_names must be an iterable of block-name strings")
        names = list(block_names)
        if not names:
            raise ValueError("block_names must not be empty")
        if not all(isinstance(name, str) for name in names):
            raise TypeError("block_names must contain only strings")
        self._cpp_obj.write_blocks(self._destination(destination), names)

    def write_parameters(
        self,
        destination: PathLike,
        parameter_ids: Iterable[ParamId],
    ) -> None:
        """Export selected parameters addressed by block and LHA identifier.

        Args:
            destination: Output filename.
            parameter_ids: Non-empty iterable of :class:`ParamId` objects.

        Raises:
            TypeError: If an entry is not a ``ParamId``.
            ValueError: If no parameter is supplied, the suffix is unsupported, or
                a requested parameter does not exist.
            RuntimeError: If HyperIso is not initialized or writing fails.
        """
        ids = list(parameter_ids)
        if not ids:
            raise ValueError("parameter_ids must not be empty")
        if not all(isinstance(parameter_id, ParamId) for parameter_id in ids):
            raise TypeError("parameter_ids must contain only ParamId instances")

        self._cpp_obj.write_parameters(
            self._destination(destination),
            [parameter_id.to_cpp() for parameter_id in ids],
        )


__all__ = ["DatabaseWriter"]
