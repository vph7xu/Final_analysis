// -----------------------------------------------------------------------------
// RunQuality.h  —  header-only helper for run‑quality flags
// -----------------------------------------------------------------------------
// * Supports two interchangeable formats:
//      1) CSV  "run , helicity_good , moller_good"
//      2) JSON {
//             "5053": {"helicity": true,  "moller": false},
//             "5054": {"helicity": true,  "moller": true }
//         }
//
// * **Missing runs are now treated as BAD** (both helicityOK & mollerOK → false)
//   so you must list every accepted run in the table.
// -----------------------------------------------------------------------------
#ifndef RUN_QUALITY_H
#define RUN_QUALITY_H

#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "json.hpp"          // nlohmann single-header JSON

class RunQuality {
public:
    // ---------------------------------------------------------
    // CSV loader: any non‑zero field counts as "good" (true)
    // ---------------------------------------------------------
    bool loadCSV(const std::string& file)
    {
        std::ifstream fin(file);
        if (!fin) return false;

        std::string line;
        while (std::getline(fin, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::istringstream iss(line);

            int run; int hel = 0; int mol = 0;
            if (!(iss >> run >> hel >> mol)) continue;
            table_[run] = { hel != 0, mol != 0 };
        }
        return true;
    }

    // ---------------------------------------------------------
    // JSON loader: expects the object form shown above.
    // Booleans OR 0/1 numbers are accepted.
    // ---------------------------------------------------------
    bool loadJSON(const std::string& file, bool ignore_comments = true)
    {
        std::ifstream f(file);
        if (!f) return false;

        const std::string txt((std::istreambuf_iterator<char>(f)),
                               std::istreambuf_iterator<char>());
        try {
            auto j = nlohmann::json::parse(txt, nullptr, true, ignore_comments);
            for (auto& [runStr, flags] : j.items()) {
                int run = std::stoi(runStr);
                bool hel = flags.value("helicity", false);
                bool mol = flags.value("moller",   false);
                table_[run] = { hel, mol };
            }
        } catch (const nlohmann::json::exception&) {
            return false;
        }
        return true;
    }

    // ---------------------------------------------------------
    // Query helpers: missing → bad (false)
    // ---------------------------------------------------------
    bool helicityOK(int run) const {
        auto it = table_.find(run);
        return it != table_.end() && it->second.helicity;
    }
    bool mollerOK(int run) const {
        auto it = table_.find(run);
        return it != table_.end() && it->second.moller;
    }

private:
    struct Flags { bool helicity = false; bool moller = false; };
    std::unordered_map<int, Flags> table_;
};

#endif // RUN_QUALITY_H

