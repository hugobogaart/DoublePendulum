//
// Created by Hugo Bogaart on 30/09/2024.
//

#ifndef SMATH_H
#define SMATH_H

#include <array>
#include <cstddef>

namespace Math
{
template <typename T, size_t N>
struct Vector : public std::array<T, N> {
        auto operator+= (const Vector<T, N>& rhs) -> Vector<T, N> &
        {
                for (size_t i = 0; i < N; i++)
                        (*this)[i] += rhs[i];
                return *this;
        }
        auto operator-= (const Vector<T, N>& rhs) -> Vector<T, N> &
        {
                for (size_t i = 0; i < N; i++)
                        (*this)[i] -= rhs[i];
                return *this;
        }
        auto operator*= (const T& rhs) -> Vector<T, N> &
        {
                for (size_t i = 0; i < N; i++)
                        (*this)[i] *= rhs;
                return *this;
        }
        auto operator/= (const T &rhs) -> Vector<T, N> &
        {
                for (size_t i = 0; i < N; i++)
                        (*this)[i] /= rhs;
                return *this;
        }
};

template <typename T, size_t N>
auto operator+(const Vector<T, N> &lhs, const Vector<T, N> &rhs) -> Vector<T, N>
{
        auto copy = lhs;
        return copy += rhs;
}
template <typename T, size_t N>
auto operator-(const Vector<T, N> &lhs, const Vector<T, N> &rhs) -> Vector<T, N>
{
        auto copy = lhs;
        return copy -= rhs;
}
template <typename T, size_t N>
auto operator*(const Vector<T, N> &lhs, const T &rhs) -> Vector<T, N>
{
        auto copy = lhs;
        return copy *= rhs;
}
template <typename T, size_t N>
auto operator*(const T &lhs, const Vector<T, N> &rhs) -> Vector<T, N>
{
        auto copy = rhs;
        return copy *= lhs;
}

template <typename T, size_t N>
auto operator/(const Vector<T, N> &lhs, const T &rhs) -> Vector<T, N>
{
        auto copy = lhs;
        return copy /= rhs;
}

} // namespace Math
#endif //SMATH_H
