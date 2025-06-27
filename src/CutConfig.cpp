#include "CutConfig.h"
#include <fstream>
#include <nlohmann/json.hpp>   // single-header JSON (header-only → no linking)

using json = nlohmann::json;

bool CutConfig::load(const std::string& filename)
{
    std::ifstream f(filename);
    if (!f.is_open()) return false;

    json j;      ///< supports comments if you use JSON5 or nlohmann v3.11+
    f >> j;
    for (auto& el : j.items())
        values_[el.key()] = el.value().get<double>();
    return true;
}

double CutConfig::operator[](const std::string& key) const
{
    auto it = values_.find(key);
    if (it == values_.end())
        throw std::runtime_error("CutConfig: key '" + key + "' not found");
    return it->second;
}

bool CutConfig::contains(const std::string& key) const
{
    return values_.count(key);
}

