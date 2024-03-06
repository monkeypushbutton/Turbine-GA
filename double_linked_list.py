#!/usr/bin/env python3

"""
    double_linked_list.py
"""

__author__    = 'Nils Olofsson'
__email__     = 'me@nilsolovsson.se'
__copyright__ = 'Copyright 2021, AllSystemsPhenomenal'
__license__   = 'MIT'

# ==============================================================================
#
# Double linked list
#
# ==============================================================================

class Node:
    """
        Node element in a DoubleLinkedList.
        Each node in a valid list is associated with a value/data element and
        with its left and right neighbor.
        [Prev. node]<--[Node]-->[Next node]
                         |
                       [Data]
    """

    def __init__(self, data):
        self.data = data
        self.prev = None
        self.next = None

class DoubleLinkedList:
    """
        A double linked list. Each element keeps a reference to both left and
        right neighbor. This allows e.g. for easy removal of elements.
        The list is circular and is usually considered traversed when the next element
        is the same element as when when we started.
    """

    def __init__(self):
        self.first = None
        self.size  = 0

    def __str__(self):
        if self.first==None:
            return '[]'
        msg = '['
        msg += str(self.first.data)
        node = self.first.next
        while node != self.first:
            msg += ', ' + str(node.data)
            node = node.next
        msg += ']'
        return msg

    def append(self, data):
        self.size += 1
        if self.first == None:
            self.first = Node(data)
            self.first.prev = self.first
            self.first.next = self.first
            return
        node = Node(data)
        last = self.first.prev
        node.prev = last
        node.next = self.first
        last.next = node
        self.first.prev = node

    def remove(self, item):
        if self.first==None:
            return
        rmv = None
        node = self.first
        if node.data == item:
            rmv = node
        node = node.next
        while not rmv and node != self.first:
            if node.data == item:
                rmv = node
            node = node.next
        if rmv:
            nxt = rmv.next
            prv = rmv.prev
            prv.next = nxt
            nxt.prev = prv
            self.size -= 1
            if rmv == self.first:
                self.first = nxt
            if rmv == self.first:
                self.first = None
        return

    def count(self):
        if self.first==None:
            return 0
        i = 1
        node = self.first.next
        while node != self.first:
            i+=1
            node = node.next
        return i

    def flatten(self):
        if self.first==None:
            return []
        l = []
        node = self.first
        l.append(node.data)
        node = self.first.next
        while node != self.first:
            l.append(node.data)
            node = node.next
        return l