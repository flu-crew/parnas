# -*- coding: utf-8 -*-
from dendropy import TaxonNamespace, Tree


class InvalidArgumentError(Exception):
    def __init__(self, name: str, value: str, message=''):
        super(InvalidArgumentError, self).__init__('Invalid argument %s: %s; %s' % (name, value, message))
        self.name = name
        self.value = value


class TreeIndexer:
    """
    Adds 'index' field to all nodes and taxa on the passed trees.
    It ensures that the taxa are indexed consistently across trees.
    """

    def __init__(self, taxon_namespace: TaxonNamespace):
        self.label_mapping = {}
        index = 0
        for taxon in taxon_namespace:
            self.label_mapping[taxon.label] = index
            index += 1

    def index_tree(self, tree: Tree):
        for leaf_node in tree.leaf_nodes():
            if leaf_node.taxon.label in self.label_mapping:
                leaf_node.taxon.index = self.label_mapping[leaf_node.taxon.label]
            else:
                print(leaf_node.taxon.label)
                print(self.label_mapping)
                raise InvalidArgumentError('tree', '', 'Input tree should be over the initially specified taxon set')

        node_id = 0
        for leaf in tree.leaf_node_iter():
            leaf.index = node_id
            node_id += 1
        for node in tree.preorder_internal_node_iter():
            node.index = node_id
            # node.annotations.add_new('id', node_id)
            node_id += 1
