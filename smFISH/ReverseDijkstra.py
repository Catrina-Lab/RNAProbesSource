from __future__ import annotations
from collections.abc import Callable, Collection, Generator


#reverse dijkstra algorithm for smFISH
class ReverseDijkstraItem:
    def __init__(self, outer_instance, index):
        self.outer_instance = outer_instance
        self.index = index
        self.max_value = 0
        self.best_next_node = None
        self.calculate_path()

    def calculate_path(self):
        self.max_value, self.best_next_node = self.outer_instance.calculate_path(self.index)

def get_path(item: ReverseDijkstraItem):
    while item is not None:
        yield item.outer_instance.get_seq_elem(item.index)
        item = item.best_next_node

class ReverseDijkstra:

    def __init__(self, sequence: Collection[any], value_mapper: Callable[any, float], can_have_path: Callable[any, any, any, bool],
                 indexer : Callable[any, int, any] = lambda seq, i: seq[i]):
        """
        Run's a reverse dijskra algorithm on a sequence (optimized to work well for smFISH) to determine the best path.
        :param sequence: the sequence to run the algorithm on. Must adhere to this rule: if for an index i and values, j>=0, k>j, l>k,:
        i+j fulfills can_have_path and i+k does not, any i+l will also not fulfill it
        :param value_mapper: from an element, get a float describing the cost increase if this item should be added. This must be never be decreasing (will always return >= 0)
        :param can_have_path: (using values a,b,c) returns if there can be a path from node a->c given b is the first path found.
            Note:
                - Should return possible starting nodes if -1 is passed as a
                - b can be none if trying to find the first match
        """
        self.sequence = sequence
        self.value_mapper = value_mapper
        self.can_have_path = can_have_path
        self.cache : list = [None] * len(sequence)
        self.indexer = indexer

    def run(self) -> tuple[float, Generator]:
        max_value, node = self.calculate_path(-1)
        return max_value, get_path(node)

    def get_seq_elem(self, index):
        return self.indexer(self.sequence, index) if index != -1 else None

    def calculate_path(self, index: int) -> tuple[int, ReverseDijkstraItem]:
        seq_item = self.get_seq_elem(index)
        i = index + 1
        best_val : ReverseDijkstraItem = None
        while i < len(self.sequence) and not self.can_have_path(seq_item, None, self.get_seq_elem(i)): i+= 1 #skip invalid

        if i < len(self.sequence):
            first_item = self.get_seq_elem(i)
            while i < len(self.sequence) and self.can_have_path(seq_item, first_item, self.get_seq_elem(i)):
                self.cache[i] = self.cache[i] or ReverseDijkstraItem(self, i)
                if best_val is None or self.cache[i].max_value > best_val.max_value:
                    best_val = self.cache[i]
                i += 1

        return ((best_val.max_value if best_val is not None else 0) +
                (0 if index == -1 else self.value_mapper(self.get_seq_elem(index))),
                best_val)
