#pragma once

#include <vector>
#include <stdexcept>

using namespace std;

template <typename T>
class SymmetricMatrix
{
private:
    vector<T> storage;
    size_t size;

public:
    SymmetricMatrix(size_t size) : size(size)
    {
        storage.reserve(size * (size + 1) / 2);
    }

    T &operator()(size_t i, size_t j)
    {
        return storage[get_index(i, j)];
    }

    const T &operator()(size_t i, size_t j) const
    {
        return storage[get_index(i, j)];
    }

private:
    size_t get_index(size_t i, size_t j) const
    {
        if (i < 0 || i >= size || j < 0 || j >= size)
            throw out_of_range("Invalid index.");

        if (i > j)
            swap(i, j);

        return (i * size) - (i * (i + 1)) / 2 + j;
    }
};