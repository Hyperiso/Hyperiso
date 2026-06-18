from pyhyperiso.phyperiso.pyhyperiso import common
from pyhyperiso.core.Common.GeneralEnum import ParameterType
from dataclasses import dataclass, field
from typing import Optional, List, Union, Set

class LhaID:
    """Python wrapper around the C++ ``LhaID`` identifier.

    ``LhaID`` represents an identifier as an *ordered* sequence of integer parts.
    This allows representing both:
      - simple IDs (single part), e.g. a PDG-like code, and
      - multi-part IDs (several parts), e.g. (block, row, column).

    The canonical string representation joins the integer parts with underscores:

    - ``[511]`` -> ``"511"``
    - ``[5, 2]`` -> ``"5_2"``
    - ``[321, 1, 3]`` -> ``"321_1_3"``

    Notes:
        - Hashing and comparisons are delegated to the underlying C++ object.
        - Converting to ``int`` keeps **only the first part** (potentially losing
          information for multi-part IDs), mirroring the C++ behavior.
    """
    
    def __init__(self, *args: Union[int, str, List[int], common.LhaID]):
        """Create a ``LhaID`` from common Python representations.

        This wrapper mirrors the C++ constructors which accept:
        - a variadic list of integral sub-IDs,
        - an underscore-separated string,
        - a vector/list of integers,
        - or an existing ``LhaID`` object.

        Args:
            *args: One of the following forms:

                - ``LhaID(cpp_id)``
                  Where ``cpp_id`` is an existing ``common.LhaID`` instance.

                - ``LhaID(n)``
                  Where ``n`` is a single integer sub-ID.

                - ``LhaID("1_2_3")``
                  Underscore-separated sub-IDs.

                - ``LhaID([1, 2, 3])``
                  List of sub-IDs.

                - ``LhaID(1, 2, 3)``
                  Multiple positional integer sub-IDs.

        Raises:
            TypeError: If the provided argument types are not supported.

        Examples:
            >>> LhaID(511).to_string()
            '511'
            >>> LhaID("5_2").get_parts()
            [5, 2]
            >>> LhaID(321, 1, 3).to_string()
            '321_1_3'
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
                raise TypeError("Unsupported argument for LhaID")
        else:
            self._cpp_obj = common.LhaID(*args)

    def to_cpp(self):
        """Return the underlying bound C++ object.

        Returns:
            common.LhaID: The wrapped pybind11 C++ instance.
        """
        return self._cpp_obj
    
    def to_string(self):
        """Return the canonical underscore-joined string representation.

        Returns:
            str: Canonical string form (e.g. ``"1_2_3"``). If the identifier has
            no parts, the C++ implementation typically returns an empty string.
        """
        return self._cpp_obj.to_string()

    def get_parts(self):
        """Return the ordered list of integer parts.

        Returns:
            List[int]: The sub-IDs making up this identifier.
        """
        return self._cpp_obj.get_parts()

    def __int__(self):
        """Convert the identifier to a Python ``int``.

        This is intended for *trivial* (single-part) IDs. If the underlying ID
        contains multiple parts, only the **first** part is returned (the C++
        implementation may emit a warning).

        Returns:
            int: The first part of the identifier.

        Raises:
            ValueError: If the identifier has no parts (depends on the C++ binding).
        """
        return int(self._cpp_obj)

    def __eq__(self, other):
        """Check equality with another ``LhaID``.

        Args:
            other (LhaID): Another wrapped identifier.

        Returns:
            bool: ``True`` if both wrap equal C++ ``LhaID`` values.
        """
        return isinstance(other, LhaID) and self._cpp_obj == other._cpp_obj

    def __repr__(self):
        """Return a debug representation.

        Returns:
            str: Debug-style string such as ``"LhaID(1_2_3)"``.
        """
        return f"LhaID({self.to_string()})"

    def __hash__(self):
        """Make the object usable as a dict key / set element.

        Returns:
            int: Hash based on the identifier parts.
        """
        return hash(tuple(self.get_parts()))
    
if __name__ == "__main__":

    print("\n🔢 Testing LhaID...")
    lid1 = LhaID(32)
    lid2 = LhaID("1_2_3")
    lid3 = LhaID([4, 5, 6])

    print(f"lid1 = {lid1}, int: {int(lid1)}, parts: {lid1.get_parts()}")
    print(f"lid2 = {lid2}, parts: {lid2.get_parts()}")
    print(f"lid3 = {lid3}, parts: {lid3.get_parts()}")
    print("✅ lid2 == LhaID('1_2_3'):", lid2 == LhaID("1_2_3"))

