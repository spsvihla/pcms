/** @file tree_node_module.h
 *  @brief Header file for tree_module.c
 *  @author Sean Svihla
 */

#ifndef TREE_NODE_MODULE_H
#define TREE_NODE_MODULE_H

#include <Python.h>

typedef struct {
    PyObject_HEAD
    PyObject *child, *sibling;
    float edge_length;
} TreeNodeStruct;

#endif // TREE_NODE_MODULE_H