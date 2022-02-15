#pragma once
#include <vector>
#include <sstream>

namespace VectorUtilities
{
    /// <summary>
    /// Utility function to convert a string to a numeric type.
    /// </summary>
    /// <typeparam name="T">Numeric type. E.g: int, double, etc</typeparam>
    /// <param name="str">Input string.</param>
    /// <returns>Converted type.</returns>
    template <typename T> T ConvertTo(const std::string& str)
    {
        std::istringstream ss(str);
        T num;
        ss >> num;
        return num;
    }

    /// <summary>
    /// Check if the vector v contains the key.
    /// </summary>
    template<typename T> bool VectorContainsQ(const std::vector<T>& v, const T& key)
    {
        return std::find(v.begin(), v.end(), key) != v.end();
    }

    template<typename T> void RemoveLower(std::vector<T>& v, T n)
    {
        auto removePositions = remove_if(v.begin(), v.end(), [&](auto const& val){ return val < n; });
        v.erase(removePositions, v.end());
    }
    template<typename T> void RemoveGreater(std::vector<T>& v, T n)
    {
        auto removePositions = remove_if(v.begin(), v.end(), [&](auto const& val){ return val > n; });
        v.erase(removePositions, v.end());
    }
    template<typename T> void RemoveGreaterOrEqual(std::vector<T>& v, T n)
    {
        auto removePositions = remove_if(v.begin(), v.end(), [&](auto const& val){ return val >= n; });
        v.erase(removePositions, v.end());
    }
    template<typename T> void Sort(std::vector<T>& v)
    {
        std::sort(v.begin(), v.end());
    }
    template<typename T> int IndexOf(const std::vector<T>& v, T n)
    {
        return static_cast<int>(upper_bound(v.begin(), v.end(), n) - v.begin());
    }
    template<typename T> int NumberOfGreater(const std::vector<T>& v, T n)
    {
        return static_cast<int>(count_if(v.begin(), v.end(),[&](auto const& val){ return val > n; }));
    }
    template<typename T> int NumberOfLower(const std::vector<T>& v, T n)
    {
        return static_cast<int>(count_if(v.begin(), v.end(),[&](auto const& val){ return val < n; }));
    }
    template<typename T> int NumberOfGreaterOrEqual(const std::vector<T>& v, T n)
    {
        return static_cast<int>(count_if(v.begin(), v.end(),[&](auto const& val){ return val >= n; }));
    }
    template<typename T> int NumberOfLowerOrEqual(const std::vector<T>& v, T n)
    {
        return static_cast<int>(count_if(v.begin(), v.end(),[&](auto const& val){ return val <= n; }));
    }
    template<typename T> T Max(const std::vector<T>& v)
    {
        return *max_element(v.begin(), v.end());
    }
    template<typename T> T Min(const std::vector<T>& v)
    {
        return *min_element(v.begin(), v.end());
    }
}