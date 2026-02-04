#include <string>

inline bool ends_with(const std::string& str, const std::string& suffix) {
    if (suffix.size() > str.size()) {
        return false; // Suffix longer than string can't match
    }
    return std::equal(suffix.rbegin(), suffix.rend(), str.rbegin());
}


