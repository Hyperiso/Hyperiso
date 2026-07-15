"""Typed names for LHA and internal parameter blocks."""

from pyhyperiso.phyperiso.pyhyperiso import common
from typing import List, Union, Set


class BlockName:
    """Python wrapper around the C++ ``BlockName`` with alias semantics.

    A ``BlockName`` is not just one string: it is a *set of aliases* that all
    refer to the same logical block (useful for legacy naming conventions or
    capitalization differences).

    Important semantic detail (from C++):
        Two block names can be considered equal if they **share at least one**
        alias (depending on the C++ operator== implementation).

    Notes:
        - ``to_string()`` returns a representative string for the block. If
          multiple aliases exist, which one is returned may be unspecified.
          Prefer ``get_alias()`` when you need all aliases.
        - ``add_alias()`` and ``to_upper()`` are chainable and return ``self``.
    """

    def __init__(self, names: Union[str, List[str], Set[str], common.BlockName]):
        """Create a ``PyBlockName`` from common Python representations.

        Args:
            names: One of:
                - ``common.BlockName``: an existing bound C++ instance.
                - ``str``: a single alias (some C++ implementations may also
                  parse slash-separated aliases like ``"A/B/C"``).
                - ``list[str]`` or ``set[str]``: multiple aliases.

        Raises:
            TypeError: If ``names`` has an unsupported type.
        """
        if isinstance(names, common.BlockName):
            self._cpp_obj = names
        elif isinstance(names, str):
            self._cpp_obj = common.BlockName(names)
        elif isinstance(names, (list, set)):
            str_set = {str(name) for name in names}
            self._cpp_obj = common.BlockName(str_set)
        else:
            raise TypeError(f"Unsupported type for BlockName init: {type(names)}")

    def to_string(self) -> str:
        """Return a representative string for this block.

        Returns:
            str: A string alias for this block. If several aliases are present,
            the chosen alias may be unspecified (mirrors C++).
        """
        return self._cpp_obj.to_string()

    def get_alias(self) -> Set[str]:
        """Return the full alias set.

        Returns:
            Set[str]: All aliases associated with this block name.
        """
        return set(self._cpp_obj.get_alias())

    def has_alias(self, alias: str) -> bool:
        """Check whether an alias is registered.

        Args:
            alias (str): Alias string to test.

        Returns:
            bool: ``True`` if ``alias`` is present in the alias set.
        """
        return self._cpp_obj.has_alias(alias)

    def add_alias(self, alias: str):
        """Add a new alias to this block name.

        Args:
            alias (str): Alias to add.

        Returns:
            PyBlockName: ``self`` (chainable).
        """
        self._cpp_obj.add_alias(alias)
        return self  # chainable

    def to_upper(self):
        """Normalize aliases to upper case.

        Returns:
            PyBlockName: ``self`` (chainable).
        """
        self._cpp_obj.to_upper()
        return self  # chainable

    def __eq__(self, other):
        """Equality comparison.

        Args:
            other (Union[PyBlockName, str]): Another block name wrapper or a string.

        Returns:
            bool: Result of the underlying C++ equality semantics. When comparing
            to a string, the binding typically checks against the alias set.
        """
        if isinstance(other, BlockName):
            return self._cpp_obj == other._cpp_obj
        elif isinstance(other, str):
            return self._cpp_obj == other
        return False

    def __ne__(self, other):
        """Negated equality."""
        return not self.__eq__(other)

    def __lt__(self, other):
        """Strict ordering.

        Args:
            other (PyBlockName): Another block name wrapper.

        Returns:
            bool: Result of the underlying C++ ordering (if provided).

        Raises:
            TypeError: If ``other`` is not a ``PyBlockName``.
        """
        if isinstance(other, BlockName):
            return self._cpp_obj < other._cpp_obj
        raise TypeError(f"Cannot compare BlockName with {type(other)}")

    def __hash__(self):
        """Hash for set/dict usage."""
        return hash(self._cpp_obj)

    def __str__(self):
        """String conversion."""
        return self.to_string()

    def __repr__(self):
        """Debug representation."""
        return f"BlockName({self.get_alias()})"
