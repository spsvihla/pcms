/** @file tree_node_module.c
 *  @brief A Python C extension TreeNode class.
 *  @author Sean Svihla
 */

#include <Python.h>
#include <structmember.h>
#include "tree_node_module.h"

extern PyTypeObject TreeNode_Type; // forward declare to type check in getsetters

/*****************************************************************************
 *                           TreeNode class definition                       *
 *****************************************************************************/

static int
TreeNode_traverse(TreeNodeStruct *self, visitproc visit, void *arg)
{
    if(self->child != NULL)
    {
        Py_VISIT(self->child);
    }

    if(self->sibling != NULL)
    {
        Py_VISIT(self->sibling);
    }

    return 0;
}

static int
TreeNode_clear(TreeNodeStruct *self)
{
    Py_CLEAR(self->child);
    Py_CLEAR(self->sibling);
    return 0;
}

static void
TreeNode_dealloc(TreeNodeStruct *self)
{
    PyObject_GC_UnTrack(self);
    TreeNode_clear(self);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject*
TreeNode_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    TreeNodeStruct *self;
    self = (TreeNodeStruct *)type->tp_alloc(type, 0); // includes Py_INCREF

    if(self == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, 
                        "Failed to allocate memory for new TreeNode object");
        return NULL;
    }

    return (PyObject *)self;
}

/* A helper function for getters to return a member of TreeNode */
static inline PyObject*
get_TreeNode_member(PyObject **member)
{
    if(*member == NULL)
    {
        Py_RETURN_NONE;
    }
    Py_INCREF(*member);
    return *member;
}

static PyObject*
TreeNode_get_child(PyObject *self, PyObject *closure)
{
    TreeNodeStruct *tree = (TreeNodeStruct *)self;
    return get_TreeNode_member(&tree->child);
}

static PyObject*
TreeNode_get_sibling(PyObject *self, PyObject *closure)
{
    TreeNodeStruct *tree = (TreeNodeStruct *)self;
    return get_TreeNode_member(&tree->sibling);
}

/* helper function for setters to set a member of TreeNode */
static inline int
set_TreeNode_member(PyObject **member, PyObject *value, 
                    const char error_message)
{
    if(value == NULL || value == Py_None)
    {
        if(*member != NULL)
        {
            Py_CLEAR(*member);
        }
    }
    else
    {
        if(!PyObject_IsInstance(value, (PyObject *)&TreeNode_Type))
        {
            PyErr_SetString(PyExc_TypeError, error_message);
            return -1;
        }

        if(*member != NULL)
        {
            Py_DECREF(*member);
        }

        Py_INCREF(value);
        *member = value;
    }
    return 0;
}

static int
TreeNode_set_child(PyObject *self, PyObject *value, PyObject *closure)
{
    TreeNodeStruct *tree = (TreeNodeStruct *)self;
    return set_TreeNode_member(&tree->child, value,
                           "Expected TreeNode type or None for child");
}

static int
TreeNode_set_sibling(PyObject *self, PyObject *value, PyObject *closure)
{
    TreeNodeStruct *tree = (TreeNodeStruct *)self;
    return set_TreeNode_member(&tree->sibling, value,
                           "Expected TreeNode type or None for sibling");
}

/* helper function for TreeNode_init to type check and set defaults */
static inline int
typecheck_default(PyObject **arg, PyTypeObject *expected_type, 
                       PyObject *default_value, const char *error_message)
{
    if(*arg == NULL || *arg == Py_None)
    {
        *arg = default_value;
    }
    else
    {
        if(!PyObject_IsInstance(*arg, (PyObject *)expected_type))
        {
            PyErr_SetString(PyExc_ValueError, error_message);
            return -1;
        }
    }
    return 0;
}

static int
TreeNode_init(TreeNodeStruct *tree, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"child", "sibling", "edge_length", NULL};
    PyObject    *child = NULL, *sibling = NULL, *edge_length = NULL,
                *default_value;

    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "|OOO", kwlist, 
                                    &child, &sibling, &edge_length))
    {
        return -1;
    }

    // set child
    if(typecheck_default(&child, &TreeNode_Type, NULL,
                      "Expected TreeNode type or None for child") == -1)
    {
        return -1;
    }
    TreeNode_set_child((PyObject *)tree, child, NULL);

    // set sibling
    if(typecheck_default(&sibling, &TreeNode_Type, NULL,
                      "Expected TreeNode type or None for sibling") == -1)
    {
        return -1;
    }
    TreeNode_set_sibling((PyObject *)tree, sibling, NULL);

    // set edge_length
    if(default_value == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "Failed to create new PyFloat object");
        return -1;
    }
    if(edge_length == NULL || edge_length == Py_None)
    {
        edge_length = PyFloat_FromDouble(0.0); // default value
    }
    else
    {
        if(!PyNumber_Check(edge_length))
        {
            PyErr_SetString(PyExc_ValueError,
                            "Expected numeric or None for edge_length");
        return -1;
    }
    }
    tree->edge_length = (float)PyFloat_AsDouble(edge_length);

    return 0;
}

static PyMemberDef TreeNode_members[] = {
    {"edge_length", T_FLOAT, offsetof(TreeNodeStruct, edge_length), 0,
     "The edge length from the parent"},
    {NULL}  /* Sentinel */
};

static PyGetSetDef TreeNode_getsetters[] = {
    {"child", (getter) TreeNode_get_child, (setter) TreeNode_set_child,
     "first child node", NULL},
    {"sibling", (getter) TreeNode_get_sibling, (setter) TreeNode_set_sibling,
     "first sibling node", NULL},
    {NULL}  /* Sentinel */
};

/*****************************************************************************
 *                         TreeNode method definitions                       *
 *****************************************************************************/

/* convenience wrapper for TreeNode_init */
static inline PyObject*
_TreeNode_init(PyObject *args, PyObject *kwargs)
{
    PyObject *child;
    int       didCreateArgs, didCreateKwargs;

    // make empty args tuple if none given
    if(args == NULL)
    {
        args = PyTuple_New(0);
        if(args == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError,
                            "Failed to create new PyTuple object");
            return NULL;
        }
        didCreateArgs = 1;
    }

    // make empty dict is none given
    if(kwargs == NULL)
    {
        kwargs = PyDict_New();
        if(kwargs == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError,
                            "Failed to create new PyDict object");
            if(didCreateArgs)
            {
                Py_CLEAR(args);
            }
            return NULL;
        }
    }

    child = TreeNode_new(&TreeNode_Type, NULL, NULL);
    if(child == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "Failed to create new child TreeNode");
        if(didCreateArgs)
        {
            Py_CLEAR(args);
        }
        if(didCreateKwargs)
        {
            Py_CLEAR(kwargs);
        }
        return NULL;
    }

    if(TreeNode_init((TreeNodeStruct *)child, args, kwargs) == -1)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "Failed to initialize child TreeNode");
        if(didCreateArgs)
        {
            Py_CLEAR(args);
        }
        if(didCreateKwargs)
        {
            Py_CLEAR(kwargs);
        }
        Py_CLEAR(child);
        return NULL;
    }

    return child;
}

static PyObject*
TreeNode_add_child(PyObject *self, PyObject *args)
{
    TreeNodeStruct *tree = (TreeNodeStruct *)self;
    PyObject       *child, *arg = NULL, *TreeNode_init_args;

    if(!PyArg_ParseTuple(args, "|O", &arg))
    {
        return NULL;
    }

    if(arg == NULL || arg == Py_None)
    {
        // arg is NULL or None, so set default values for new TreeNode members   

        child = _TreeNode_init(NULL, NULL);
        if(child == NULL)
        {
            // PyErr_SetString in _TreeNode_init
            return NULL;
        }
    }
    else if(PyObject_IsInstance(arg, (PyObject *)&TreeNode_Type))
    {
        // the new child is given, so use it
        child = arg;
    }
    else if(PyNumber_Check(arg))
    {
        // arg is a numeric object, so create a PyTuple "TreeNode_init_args" to 
        // pass to TreeNode_init

        TreeNode_init_args = PyTuple_New(3);
        if(TreeNode_init_args == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError,
                            "Failed to create new PyTuple object");
            return NULL;
        }
        PyTuple_SET_ITEM(TreeNode_init_args, 0, Py_None);
        PyTuple_SET_ITEM(TreeNode_init_args, 1, Py_None);
        PyTuple_SET_ITEM(TreeNode_init_args, 2, arg);

        child = _TreeNode_init(TreeNode_init_args, NULL);
        if(child == NULL)
        {
            // PyErr_SetString in _TreeNode_init
            Py_CLEAR(TreeNode_init_args);
            return NULL;
        }
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, 
                        "Expected TreeNode object, numeric, or None");
        return NULL;
    }

    // add to tree->child if it is NULL; otherwise add to the last sibling of 
    // treee->child
    if(tree->child == NULL)
    {
        TreeNode_set_child((PyObject *)tree, child, NULL);
    }
    else
    {
        tree = (TreeNodeStruct *)tree->child;
        while(tree->sibling != NULL)
        {
            tree = (TreeNodeStruct *)tree->sibling;
        }
        TreeNode_set_sibling((PyObject *)tree, child, NULL);   
    }

    Py_CLEAR(TreeNode_init_args);
    return child;
}

static PyMethodDef TreeNode_methods[] = {
    {"add_child", TreeNode_add_child, METH_VARARGS, 
     "Add a child node to the tree"},
    {NULL}  /* Sentinel */
};

PyTypeObject TreeNode_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "tree_node.TreeNode",
    .tp_doc = "A TreeNode data structure",
    .tp_basicsize = sizeof(TreeNodeStruct),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .tp_new = TreeNode_new,
    .tp_init = (initproc) TreeNode_init,
    .tp_dealloc = (destructor) TreeNode_dealloc,
    .tp_traverse = (traverseproc) TreeNode_traverse,
    .tp_clear = (inquiry) TreeNode_clear,
    .tp_members = TreeNode_members,
    .tp_methods = TreeNode_methods,
    .tp_getset = TreeNode_getsetters,
};

/*****************************************************************************
 *                     top-level method definitions                          *
 *****************************************************************************/

static PyMethodDef TreeNodeMethods[] = {
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

/*****************************************************************************
 *                            module definition                              *
 *****************************************************************************/

static PyModuleDef treenodemodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "tree_node",
    .m_doc = "The tree_node.TreeNode class and helper functions",
    .m_size = -1,
    .m_methods = TreeNodeMethods
};

PyMODINIT_FUNC PyInit_tree_node(void)
{
    PyObject *m;

    if(PyType_Ready(&TreeNode_Type) < 0)
    {
        return NULL;
    }

    m = PyModule_Create(&treenodemodule);
    if(m == NULL)
    {
        return NULL;
    }

    Py_INCREF(&TreeNode_Type);
    if(PyModule_AddObject(m, "TreeNode", (PyObject *)&TreeNode_Type) < 0)
    {
        Py_DECREF(&TreeNode_Type);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
