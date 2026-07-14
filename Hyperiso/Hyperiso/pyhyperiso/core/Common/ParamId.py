from pyhyperiso.phyperiso.pyhyperiso import common
from pyhyperiso.core.Common.GeneralEnum import ParameterType
from dataclasses import dataclass, field
from typing import Optional, List, Union
from pyhyperiso.core.Common.BlockName import BlockName
from pyhyperiso.core.Common.LhaID import LhaID


@dataclass
class ParamId:
    """Python wrapper around the C++ ``ParamId`` composite identifier.

    A ``ParamId`` uniquely identifies a parameter in LHA/SLHA-like structures by:
      - an optional semantic parameter type (``ParameterType``),
      - a block name (``BlockName`` with alias support),
      - an index or multi-index within that block (``LhaID``).

    This mirrors the C++ struct:

    - ``type`` is optional (can be unset / ``None``),
    - ``block`` is a block identifier with aliases,
    - ``code`` is an ``LhaID`` (single-part or multi-part).

    Attributes:
        type (Optional[ParameterType]): Optional high-level parameter category.
            ``None`` means "unset" (equivalent to ``std::nullopt`` in C++).
        block (Union[str, BlockName]): Block identifier. The dataclass field is
            declared as ``str`` for convenience, but after initialization it is
            normalized to a ``BlockName`` instance.
        code (Union[LhaID, int, str, List[int]]): Parameter index / multi-index.
            After initialization it is normalized to a ``LhaID`` instance.

    Notes:
        - The default/null ``ParamId`` corresponds to:
          ``type=None``, ``block="NULL"``, ``code=0``.
        - This class keeps an underlying C++ object in ``_cpp_obj``.
    """

    type: Optional[ParameterType] = None
    block: str = "NULL"
    code: Union[LhaID, int, str, List[int]] = field(default_factory=lambda: LhaID(0))
    _cpp_obj: common.ParamId = field(init=False, repr=False)

    def __post_init__(self):
        """Build and normalize the underlying C++ ``ParamId`` after init.

        This method:
          1) Normalizes ``code`` to ``PyLhaID`` if needed.
          2) Normalizes ``block`` to ``PyBlockName`` if needed.
          3) Constructs the C++ ``ParamId`` depending on whether ``type`` is set
             and whether the (block, code) pair corresponds to the default/null
             sentinel.

        Raises:
            TypeError: If ``block`` or ``code`` cannot be converted to the expected
                wrapper types.
        """
        if not isinstance(self.code, LhaID):
            self.code = LhaID(self.code)

        if not isinstance(self.block, BlockName):
            self.block = BlockName(self.block)

        if self.type is not None:
            self._cpp_obj = common.ParamId(self.type.value, self.block._cpp_obj, self.code._cpp_obj)
        else:
            if str(self.block) == "NULL" and int(self.code) == 0:
                self._cpp_obj = common.ParamId()
            else:
                self._cpp_obj = common.ParamId(self.block._cpp_obj, self.code._cpp_obj)

        self.type = ParameterType(self._cpp_obj.type) if self._cpp_obj.type is not None else None
        self.block = BlockName(self._cpp_obj.block)
        self.code = LhaID(self._cpp_obj.code)

    def set_parameter_type(self, param_type: ParameterType):
        """Set or overwrite the parameter type.

        Args:
            param_type (ParameterType): New parameter type to assign.
        """
        self._cpp_obj.set_parameter_type(param_type.value)
        self.type = param_type

    def to_dict(self):
        """Serialize this identifier into a Python dictionary.

        The output is intended for logging/debugging and lightweight serialization.

        Returns:
            dict: A mapping with keys:
              - ``"type"``: type name as ``str`` or ``None`` if unset,
              - ``"block"``: the current ``PyBlockName`` object,
              - ``"code"``: canonical string representation of the ``LhaID``.

        Notes:
            If you need a JSON-serializable dict, you may want:
            ``{"block": str(self.block)}`` instead of returning the object.
        """
        return {
            "type": self.type.name if self.type else None,
            "block": self.block,
            "code": self.code.to_string(),
        }

    @classmethod
    def from_cpp(cls, cpp_obj: common.ParamId) -> "ParamId":
        """Wrap an existing bound C++ ``ParamId`` instance.

        Args:
            cpp_obj (common.ParamId): C++ object coming from pybind.

        Returns:
            PyParamId: A Python wrapper around the provided C++ object.
        """
        instance = cls()
        instance._cpp_obj = cpp_obj
        instance.type = ParameterType(cpp_obj.type)
        instance.block = BlockName(cpp_obj.block)
        instance.code = LhaID(cpp_obj.code)
        return instance

    def to_cpp(self):
        """Return the underlying bound C++ object.

        Returns:
            common.ParamId: The wrapped pybind11 C++ instance.
        """
        return self._cpp_obj

    def __repr__(self):
        """Return a debug representation.

        Returns:
            str: Debug-style string including type, block and code.
        """
        return f"ParamId(type={self.type}, block='{self.block}', code={self.code})"

    def __eq__(self, other):
        if not isinstance(other, ParamId):
            return NotImplemented
        return (
            self.type,
            str(self.block),
            self.code.to_string(),
        ) == (
            other.type,
            str(other.block),
            other.code.to_string(),
        )

    def __hash__(self):
        return hash(
            (
                self.type,
                str(self.block),
                self.code.to_string(),
            )
        )
