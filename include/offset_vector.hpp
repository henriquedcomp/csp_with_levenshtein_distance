#pragma once

#include <vector>

using namespace std;

template <typename T>
class OffsetVector
{
private:
    vector<T> data;
    int offset;

public:
    // [p, q]
    OffsetVector(int p, int q) : offset(p)
    {
        if (q >= p)
            data.resize(q - p + 1);
    }

    T &operator[](int index)
    {
        return data[index - offset];
    }

    const T &operator[](int index) const
    {
        return data[index - offset];
    }

    size_t size()
    {
        return data.size();
    }

    class Iterator
    {
    private:
        T *ptr;

    public:
        explicit Iterator(T *p) : ptr(p) {}

        T &operator*() { return *ptr; }

        Iterator &operator++()
        {
            ptr++;
            return *this;
        }

        bool operator!=(const Iterator &other) const
        {
            return ptr != other.ptr;
        }
    };

    Iterator begin() { return Iterator(data.data()); }
    Iterator end() { return Iterator(data.data() + data.size()); }
};