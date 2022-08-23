"""
Find connected components.
"""
from abc import abstractmethod
from collections import OrderedDict
from typing import TypeVar, Generic, Optional, Iterable
import typing

if typing.TYPE_CHECKING:
    from typing_extensions import Protocol
else:
    Protocol = object


T = TypeVar("T")
C = TypeVar("C", bound="Comparable")


class Comparable(Protocol):
    @abstractmethod
    def __lt__(self: C, other: C) -> bool:
        pass


class Node(Generic[C]):
    def __init__(self, value: C, parent: Optional["Node"]):
        self.value = value
        self.parent = parent

    def __repr__(self):
        return f"Node(value={self.value}, parent={self.parent})"


class ComponentFinder(Generic[C]):
    """
    Find connected components. A ComponentFinder is initialized with a list of
    values. These are initially partitioned such that each value is in a
    separate set. By calling merge(x, y), the two sets containing values x and
    y are merged. Calling find(x) returns a "representative" value of the set
    that value x is in. x and y are in the same set iff find(x) == find(y).
    The representative is always the minimum value of the set.

    This implements a variant of the Union-Find algorithm, but without the
    "union by rank" strategy since we want the smallest node to be the
    representative. It could perhaps be optimized, but this function is not
    the current bottleneck.
    """

    def __init__(self, values: Iterable[C]):
        self.nodes = {value: Node(value, None) for value in values}

    def merge(self, x: C, y: C) -> None:
        assert x != y
        x_root = self._find_node(x)
        y_root = self._find_node(y)

        if x_root is y_root:
            return

        # Merge while making sure that the node with the smaller value is the
        # new parent.
        if x_root.value < y_root.value:
            y_root.parent = x_root
        else:
            x_root.parent = y_root

    def _find_node(self, value: C) -> Node[C]:
        node = root = self.nodes[value]
        while root.parent is not None:
            root = root.parent

        # compression path
        while node.parent is not None:
            node.parent, node = root, node.parent
        return root

    def find(self, value: C) -> C:
        """
        Return which component x belongs to, identified by the smallest value.
        """
        return self._find_node(value).value

    def print(self):
        for x in sorted(self.nodes):
            print(x, ":", self.nodes[x], "is represented by", self._find_node(x))


class CyclicGraphError(Exception):
    pass


class Graph:
    """Directed graph that can sort topologically"""

    def __init__(self):
        # map node to a list of neighbors
        self._neighbors = OrderedDict()

    def add_edge(self, node1, node2):
        """The edge is directed from node1 to node2"""
        if node1 not in self._neighbors:
            self._neighbors[node1] = []
        self._neighbors[node1].append(node2)
        if node2 not in self._neighbors:
            self._neighbors[node2] = []

    def toposorted(self):
        """
        Return nodes of the graph sorted topologically.
        For all edges u -> v that the graph has, node v will appear
        before node u.
        """
        order = []
        colors = {node: "white" for node in self._neighbors}

        def visit(node):
            assert colors[node] == "white"
            colors[node] = "gray"
            for neighbor in self._neighbors[node]:
                if colors[neighbor] == "white":
                    visit(neighbor)
                elif colors[neighbor] == "gray":
                    raise CyclicGraphError(f"Cycle involving {node!r} and {neighbor!r} detected")
            order.append(node)
            colors[node] = "black"

        for node in self._neighbors:
            if colors[node] == "white":
                visit(node)
        return order
