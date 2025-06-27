#ifndef CUT_CONFIG_H
#define CUT_CONFIG_H

#include <string>
#include <unordered_map>

/** Loads numeric cut thresholds from a text/JSON/YAML file */
class CutConfig {
public:
    bool load(const std::string& filename);        ///< returns false if file not found
    double operator[](const std::string& key) const; ///< convenient accessor
    bool   contains(const std::string& key) const;

private:
    std::unordered_map<std::string, double> values_;
};

#endif  // CUT_CONFIG_H

