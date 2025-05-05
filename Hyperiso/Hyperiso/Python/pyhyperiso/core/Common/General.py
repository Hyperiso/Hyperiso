from pyhyperiso.phyperiso.pyhyperiso import common
from pyhyperiso.core.Common.GeneralEnum import ParameterType
from dataclasses import dataclass, field
from typing import Optional, List, Union

class PyLhaID:
    """Python wrapper for the C++ LhaID class.

    This class allows flexible construction and usage of a LhaID object
    using Python-native types (int, str, list).
    """
    def __init__(self, *args: Union[int, str, List[int], common.LhaID]):
        """Initializes a PyLhaID object.

        Args:
            *args: Can be one of:
                - int: a single sub-ID.
                - str: a string of sub-IDs separated by underscores (e.g., "1_2_3").
                - list[int]: a list of sub-IDs.
                - common.LhaID: a pre-existing wrapped LhaID object.
        """
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, common.LhaID):
                self._cpp_obj = arg
            elif isinstance(arg, int):
                self._cpp_obj = common.LhaID(arg)
            elif isinstance(arg, str):
                self._cpp_obj = common.LhaID(arg)
            elif isinstance(arg, list):
                self._cpp_obj = common.LhaID(arg)
            else:
                raise TypeError("Unsupported argument for PyLhaID")
        else:
            self._cpp_obj = common.LhaID(*args)

    def to_string(self):
        """Returns a string representation of the LhaID.

        Returns:
            str: The string format (e.g., "1_2_3").
        """
        return self._cpp_obj.to_string()

    def get_parts(self):
        """Gets the list of sub-identifiers (parts).

        Returns:
            List[int]: List of sub-IDs.
        """
        return self._cpp_obj.get_parts()

    def __int__(self):
        """Casts the LhaID to an integer (only uses the first sub-ID).

        Returns:
            int: The first part of the LhaID.

        Raises:
            ValueError: If the LhaID has no parts.
        """
        return int(self._cpp_obj)

    def __eq__(self, other):
        """Checks equality with another PyLhaID.

        Args:
            other (PyLhaID): The object to compare.

        Returns:
            bool: True if equal.
        """
        return isinstance(other, PyLhaID) and self._cpp_obj == other._cpp_obj

    def __repr__(self):
        """Returns a string representation for debugging.

        Returns:
            str: Debug-style string.
        """
        return f"PyLhaID({self.to_string()})"


@dataclass
class PyParamId:
    """Python wrapper for the C++ ParamId class.

    Encapsulates the ParamId structure with Python-friendly access and construction.
    """
    type: Optional[ParameterType] = None
    block: str = "NULL"
    code: Union[PyLhaID, int, str, List[int]] = field(default_factory=lambda: PyLhaID(0))
    _cpp_obj: common.ParamId = field(init=False, repr=False)

    def __post_init__(self):
        """Initializes the underlying C++ ParamId object after dataclass init."""
        if not isinstance(self.code, PyLhaID):
            self.code = PyLhaID(self.code)

        cpp_type = self.type.value if self.type is not None else None

        if cpp_type is not None:
            self._cpp_obj = common.ParamId(cpp_type, self.block, self.code._cpp_obj)
        else:
            if self.block == "NULL" and int(self.code) == 0:
                self._cpp_obj = common.ParamId()
            else:
                self._cpp_obj = common.ParamId(self.block, self.code._cpp_obj)

        cpp_value = self._cpp_obj.type
        self.type = ParameterType(cpp_value) if cpp_value is not None else None
        self.block = self._cpp_obj.block
        self.code = PyLhaID(self._cpp_obj.code)

    def set_parameter_type(self, param_type: ParameterType):
        """Sets the parameter type in the underlying C++ object.

        Args:
            param_type (ParameterType): The parameter type to set.
        """
        self._cpp_obj.set_parameter_type(param_type.value)
        self.type = param_type

    def to_dict(self):
        """Serializes the ParamId to a dictionary.

        Returns:
            dict: A dict with keys 'type', 'block', and 'code'.
        """
        return {
            "type": self.type.name if self.type else None,
            "block": self.block,
            "code": self.code.to_string()
        }

    def __repr__(self):
        """Returns a string representation for debugging.

        Returns:
            str: Debug-style string.
        """
        return f"PyParamId(type={self.type}, block='{self.block}', code={self.code})"