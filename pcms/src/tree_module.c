/** @file tree_module.c
 *  @brief A Python C extension Tree class.
 *  @author Sean Svihla
 */

#include <Python.h>
#include <structmember.h>
#include "tree_module.h"

extern PyTypeObject Tree_Type; // forward declare to type check in getsetters

/*****************************************************************************
 *                           Tree class definition                           *
 *****************************************************************************/

static int
Tree_traverse(TreeStruct *self, visitproc visit, void *arg)
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
Tree_clear(TreeStruct *self)
{
    Py_CLEAR(self->child);
    Py_CLEAR(self->sibling);
    return 0;
}

static void
Tree_dealloc(TreeStruct *self)
{
    PyObject_GC_UnTrack(self);
    Tree_clear(self);
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject*
Tree_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    TreeStruct *self;
    self = (TreeStruct *)type->tp_alloc(type, 0); // includes Py_INCREF

    if(self == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, 
                        "Failed to allocate memory for new Tree object");
        return NULL;
    }

    return (PyObject *)self;
}

/* A helper function for getters to return a member of Tree */
static inline PyObject*
get_tree_member(PyObject **member)
{
    if(*member == NULL)
    {
        Py_RETURN_NONE;
    }
    Py_INCREF(*member);
    return *member;
}

static PyObject*
Tree_get_child(PyObject *self, PyObject *closure)
{
    TreeStruct *tree = (TreeStruct *)self;
    return get_tree_member(&tree->child);
}

static PyObject*
Tree_get_sibling(PyObject *self, PyObject *closure)
{
    TreeStruct *tree = (TreeStruct *)self;
    return get_tree_member(&tree->sibling);
}

/* A helper function for setters to set a member of Tree */
static inline int
set_tree_member(PyObject **member, PyObject *value, const char error_message)
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
        if(!PyObject_IsInstance(value, (PyObject *)&Tree_Type))
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
Tree_set_child(PyObject *self, PyObject *value, PyObject *closure)
{
    TreeStruct *tree = (TreeStruct *)self;
    return set_tree_member(&tree->child, value,
                           "Expected Tree type or None for child");
}

static int
Tree_set_sibling(PyObject *self, PyObject *value, PyObject *closure)
{
    TreeStruct *tree = (TreeStruct *)self;
    return set_tree_member(&tree->sibling, value,
                           "Expected Tree type or None for sibling");
}

/* A helper function for Tree_init to check arguments and set default values */
static inline int
check_arg_with_default(PyObject **arg, PyTypeObject *expected_type, 
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
Tree_init(TreeStruct *tree, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"child", "sibling", "edge_length", NULL};
    PyObject    *child = NULL, *sibling = NULL, *edge_length = NULL,
                *default_value;

    if(!PyArg_ParseTupleAndKeywords(args, kwargs, "|OOO", kwlist, 
                                    &child, &sibling, &edge_length))
    {
        return -1;
    }

    if(check_arg_with_default(&child, &Tree_Type, NULL,
                      "Expected Tree type or None for child") == -1)
    {
        return -1;
    }
    Tree_set_child((PyObject *)tree, child, NULL);

    if(check_arg_with_default(&sibling, &Tree_Type, NULL,
                      "Expected Tree type or None for sibling") == -1)
    {
        return -1;
    }
    Tree_set_sibling((PyObject *)tree, sibling, NULL);

    default_value = PyFloat_FromDouble(0.0);
    if(default_value == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "Failed to create new PyFloat object");
        return -1;
    }
    if(check_arg_with_default(&edge_length, &PyFloat_Type, default_value,
                      "Expected float or None for edge_length") == -1)
    {
        return -1;
    }
    tree->edge_length = (float)PyFloat_AsDouble(edge_length);

    return 0;
}

static PyMemberDef Tree_members[] = {
    {"edge_length", T_FLOAT, offsetof(TreeStruct, edge_length), 0,
     "The edge length from the parent"},
    {NULL}  /* Sentinel */
};

static PyGetSetDef Tree_getsetters[] = {
    {"child", (getter) Tree_get_child, (setter) Tree_set_child,
     "first child node", NULL},
    {"sibling", (getter) Tree_get_sibling, (setter) Tree_set_sibling,
     "first sibling node", NULL},
    {NULL}  /* Sentinel */
};

/*****************************************************************************
 *                         Tree method definitions                           *
 *****************************************************************************/

/* A helper function for Tree_add_child to create new children */
static inline PyObject*
create_child(PyObject *args, PyObject *kwargs)
{
    PyObject *child;

    child = Tree_new(&Tree_Type, NULL, NULL);
    if(child == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "Failed to create new child tree");
        return NULL;
    }

    if(Tree_init((TreeStruct *)child, args, kwargs) == -1)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "Failed to initialize child tree");
        Py_CLEAR(child);
        return NULL;
    }

    return child;
}

/**
 * @brief Add a child the Tree
 * 
 * The function adds a child to self. If an instance of Tree is given, then
 * that object is added as the child; otherwise, a new one is created. If
 * an int is given, then it is assumed to be the edge length of the newly 
 * created child. If no arguments are given, a new child is created with the
 * default edge length. Only one argument may be given at a time.
 * 
 * If self->child is NULL, then the child is added to self->child; otherwise,
 * the child is added as a sibling to self->child by following self->child
 * and then recursively child->sibling until a NULL sibling is reached.
 * 
 * @param self An instance of the Tree class
 * @param args PyTuple of positional arguments
 * @return The newly added (or created) child 
 */
static PyObject*
Tree_add_child(PyObject *self, PyObject *args)
{
    TreeStruct *tree = (TreeStruct *)self;
    PyObject   *child, *arg = NULL, *Tree_init_args, *Tree_init_kwargs;

    if(!PyArg_ParseTuple(args, "|O", &arg))
    {
        return NULL;
    }

    Tree_init_kwargs = PyDict_New();
    if(Tree_init_kwargs == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError,
                        "Failed to create new PyTuple object");
        return NULL;
    }

    if(arg == NULL || arg == Py_None)
    {
        // arg is NULL or None, so create a empty args tuple so that the new
        // child has default values for its members   
         
        Tree_init_args = PyTuple_New(0);
        if(Tree_init_args == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError,
                            "Failed to create new PyTuple object");
            Py_CLEAR(Tree_init_kwargs);
            return NULL;
        }

        child = create_child(Tree_init_args, Tree_init_kwargs);
        if(child == NULL)
        {
            // PyErr_SetString in create_child
            Py_CLEAR(Tree_init_args);
            Py_CLEAR(Tree_init_kwargs);
            return NULL;
        }
    }
    else if(PyObject_IsInstance(arg, (PyObject *)&Tree_Type))
    {
        // the new child is given, so use it
        child = arg;
    }
    else if(PyFloat_Check(arg))
    {
        // arg is a PyFloat object, so create a PyTuple "Tree_init_args" to 
        // pass to Tree_init

        Tree_init_args = PyTuple_New(3);
        if(Tree_init_args == NULL)
        {
            PyErr_SetString(PyExc_RuntimeError,
                            "Failed to create new PyTuple object");
            Py_CLEAR(Tree_init_kwargs);
            return NULL;
        }
        PyTuple_SET_ITEM(Tree_init_args, 0, Py_None);
        PyTuple_SET_ITEM(Tree_init_args, 1, Py_None);
        PyTuple_SET_ITEM(Tree_init_args, 2, arg);

        child = create_child(Tree_init_args, Tree_init_kwargs);
        if(child == NULL)
        {
            // PyErr_SetString in create_child
            Py_CLEAR(Tree_init_args);
            Py_CLEAR(Tree_init_kwargs);
            return NULL;
        }
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, 
                        "Expected Tree object, float, or None");
        Py_CLEAR(Tree_init_kwargs);
        return NULL;
    }

    // add to tree->child if it is NULL; otherwise add to the last sibling of 
    // treee->child
    if(tree->child == NULL)
    {
        Tree_set_child((PyObject *)tree, child, NULL);
    }
    else
    {
        tree = (TreeStruct *)tree->child;
        while(tree->sibling != NULL)
        {
            tree = (TreeStruct *)tree->sibling;
        }
        Tree_set_sibling((PyObject *)tree, child, NULL);   
    }

    Py_CLEAR(Tree_init_args);
    Py_CLEAR(Tree_init_kwargs);
    return child;
}

static PyMethodDef Tree_methods[] = {
    {"add_child", Tree_add_child, METH_VARARGS, "Add a child node to the tree"},
    {NULL}  /* Sentinel */
};

PyTypeObject Tree_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "tree.Tree",
    .tp_doc = "A Tree data structure",
    .tp_basicsize = sizeof(TreeStruct),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,
    .tp_new = Tree_new,
    .tp_init = (initproc) Tree_init,
    .tp_dealloc = (destructor) Tree_dealloc,
    .tp_traverse = (traverseproc) Tree_traverse,
    .tp_clear = (inquiry) Tree_clear,
    .tp_members = Tree_members,
    .tp_methods = Tree_methods,
    .tp_getset = Tree_getsetters,
};

/*****************************************************************************
 *                     top-level method definitions                          *
 *****************************************************************************/

static PyMethodDef TreeMethods[] = {
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

/*****************************************************************************
 *                            module definition                              *
 ******************************************************************************/

static PyModuleDef treemodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "tree",
    .m_doc = "The tree.Tree class and helper functions",
    .m_size = -1,
    .m_methods = TreeMethods
};

PyMODINIT_FUNC PyInit_tree(void)
{
    PyObject *m;

    if(PyType_Ready(&Tree_Type) < 0)
    {
        return NULL;
    }

    m = PyModule_Create(&treemodule);
    if(m == NULL)
    {
        return NULL;
    }

    Py_INCREF(&Tree_Type);
    if(PyModule_AddObject(m, "Tree", (PyObject *)&Tree_Type) < 0)
    {
        Py_DECREF(&Tree_Type);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
