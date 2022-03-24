#pragma once
#include <vector>
#include <sstream>
#include <numeric>

namespace VectorUtilities
{
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
    template<typename T> T Total(std::vector<T>& v)
    {
        return std::accumulate(v.begin(), v.end(), T());
    }
    template<typename T> void Insert(std::vector<T>& v1, const std::vector<T>& v2)
    {
        v1.insert(v1.end(), v2.begin(), v2.end());
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
    template<typename T> int NumberInInterval(const std::vector<T>& v, T min, T max)
    {
        return static_cast<int>(count_if(v.begin(), v.end(),[&](auto const& val){ return (val >= min && val <= max); }));
    }
    template<typename T> T Max(const std::vector<T>& v)
    {
        return *max_element(v.begin(), v.end());
    }
    template<typename T> int IndexOfMax(const std::vector<T>& v)
    {
        return static_cast<int>(max_element(v.begin(), v.end()) - v.begin());
    }
    template<typename T> T Min(const std::vector<T>& v)
    {
        return *min_element(v.begin(), v.end());
    }
    template<typename T> int IndexOfMin(const std::vector<T>& v)
    {
        return static_cast<int>(min_element(v.begin(), v.end()) - v.begin());
    }
}