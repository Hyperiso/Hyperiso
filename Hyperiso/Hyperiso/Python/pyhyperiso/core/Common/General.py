from pyhyperiso.phyperiso.pyhyperiso import common
from pyhyperiso.core.Common.GeneralEnum import ParameterType
from dataclasses import dataclass, field
from typing import Optional, List, Union, Set

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

class PyBlockName:
    """Python wrapper for the C++ BlockName class.

    Accepts various input types (string, list, set, C++ BlockName).
    """
    def __init__(self, names: Union[str, List[str], Set[str], common.BlockName]):
        if isinstance(names, common.BlockName):
            self._cpp_obj = names
        elif isinstance(names, str):
            self._cpp_obj = common.BlockName(names)
        elif isinstance(names, (list, set)):
            str_set = {str(name) for name in names}
            self._cpp_obj = common.BlockName(str_set)
        else:
            raise TypeError(f"Unsupported type for PyBlockName init: {type(names)}")

    def to_string(self) -> str:
        return self._cpp_obj.to_string()

    def get_alias(self) -> Set[str]:
        return set(self._cpp_obj.get_alias())

    def has_alias(self, alias: str) -> bool:
        return self._cpp_obj.has_alias(alias)

    def add_alias(self, alias: str):
        self._cpp_obj.add_alias(alias)
        return self  # chainable

    def to_upper(self):
        self._cpp_obj.to_upper()
        return self  # chainable

    def __eq__(self, other):
        if isinstance(other, PyBlockName):
            return self._cpp_obj == other._cpp_obj
        elif isinstance(other, str):
            return self._cpp_obj == other
        return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, PyBlockName):
            return self._cpp_obj < other._cpp_obj
        raise TypeError(f"Cannot compare PyBlockName with {type(other)}")

    def __hash__(self):
        return hash(self._cpp_obj)

    def __str__(self):
        return self.to_string()

    def __repr__(self):
        return f"PyBlockName({self.get_alias()})"
    

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
        
        if not isinstance(self.block, PyBlockName):
            self.block = PyBlockName(self.block)

        if self.type is not None:
            self._cpp_obj = common.ParamId(self.type.value, self.block._cpp_obj, self.code._cpp_obj)
        else:
            if str(self.block) == "NULL" and int(self.code) == 0:
                self._cpp_obj = common.ParamId()
            else:
                self._cpp_obj = common.ParamId(self.block._cpp_obj, self.code._cpp_obj)

        self.type = ParameterType(self._cpp_obj.type) if self._cpp_obj.type is not None else None
        self.block = PyBlockName(self._cpp_obj.block)
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

    @classmethod
    def from_cpp(cls, cpp_obj: common.ParamId) -> "PyParamId":
        """Constructs a PyParamId from a native C++ ParamId."""
        instance = cls()
        instance._cpp_obj = cpp_obj
        instance.type = ParameterType(cpp_obj.type)
        instance.block = PyBlockName(cpp_obj.block)
        instance.code = PyLhaID(cpp_obj.code)
        return instance
    
    def __repr__(self):
        """Returns a string representation for debugging.

        Returns:
            str: Debug-style string.
        """
        return f"PyParamId(type={self.type}, block='{self.block}', code={self.code})"
    
    
if __name__ == "__main__":
    print("🔧 Testing PyBlockName...")
    b1 = PyBlockName("MASS")
    b2 = PyBlockName(["mass", "MASS"])
    b3 = PyBlockName({"mass", "other"})

    print(f"b1: {b1}")
    print(f"b2: {b2}")
    print(f"b3 (before to_upper): {b3}")

    print("✅ b1 == 'MASS' :", b1 == "MASS")
    print("✅ b2 == b1     :", b2 == b1)
    print("✅ 'mass' in b3.get_alias():", "mass" in b3.get_alias())

    b3.to_upper()
    print(f"b3 (after to_upper): {b3}")

    b3.add_alias("EXTRA")
    print(f"b3 (after add_alias): {b3}")
    print("✅ has_alias('EXTRA'):", b3.has_alias("EXTRA"))

    print("\n🔢 Testing PyLhaID...")
    lid1 = PyLhaID(32)
    lid2 = PyLhaID("1_2_3")
    lid3 = PyLhaID([4, 5, 6])

    print(f"lid1 = {lid1}, int: {int(lid1)}, parts: {lid1.get_parts()}")
    print(f"lid2 = {lid2}, parts: {lid2.get_parts()}")
    print(f"lid3 = {lid3}, parts: {lid3.get_parts()}")
    print("✅ lid2 == PyLhaID('1_2_3'):", lid2 == PyLhaID("1_2_3"))

    print("\n🧬 Testing PyParamId...")
    pid1 = PyParamId(ParameterType.SM, "MASS", 32)
    pid2 = PyParamId(block="YUKAWA", code="3_3")  # sans type
    pid3 = PyParamId(type=None, block=PyBlockName("QNUMBERS"), code=PyLhaID([6, 2]))

    print(f"pid1 = {pid1}")
    print(f"pid2 = {pid2}")
    print(f"pid3 = {pid3}")

    print("\n📤 Serialized dicts:")
    print("pid1.to_dict() =", pid1.to_dict())
    print("pid2.to_dict() =", pid2.to_dict())
    print("pid3.to_dict() =", pid3.to_dict())

    print("\n🛠️ Testing set_parameter_type...")
    pid2.set_parameter_type(ParameterType.WILSON)
    print(f"pid2 (after setting type to WILSON) = {pid2}")

    print("\n✅ Equality check:")
    print("pid1 == pid1:", pid1 == pid1)
    print("pid1 != pid2:", pid1 != pid2)