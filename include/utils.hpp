#ifndef UTILS_HPP
#define UTILS_HPP

// standard library includes
#include <memory>
#include <string>
#include <vector>

// Pybind11 includes
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;


template <typename T>
inline py::array_t<T> 
std_vec2py_array_t(const std::vector<T>& arr) 
{
    if(arr.empty()) 
    {
        return py::array_t<T>();
    }
    py::array_t<T> py_arr(arr.size());
    std::memcpy(py_arr.mutable_data(), arr.data(), arr.size() * sizeof(T));
    return py_arr;
}

template <typename T>
inline py::list 
std_vecvec2py_list_array_t(const std::vector<std::vector<T>>& vv)
{
    py::list result;
    for(const auto& vec : vv)
    {
        py::array_t<T> arr(vec.size(), vec.data());
        result.append(arr);
    }
    return result;
}

#endif // UTILS_HPP