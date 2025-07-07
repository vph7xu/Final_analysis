// -----------------------------------------------------------------------------
// PhysicsAsymmetryMulti.h  —  combine raw asym with *N* correction files
// -----------------------------------------------------------------------------
// Motivation: previous PhysicsAsymmetry assumed the three ingredients lived in
// single files.  In practice each correction (acceptance, pion, inelastic …)
// may be stored in its own ROOT or TXT file.  This class lets you register an
// arbitrary list of files, each containing a TNamed("<key>","value,error") or
// the simple two‑number TXT fallback.
//
// Usage:
//
//     PhysicsAsymmetryMulti phys;
//     phys.addSource("raw_asym.root" , "raw_asym");
//     phys.addSource("Aacc.root"     , "Aacc");
//     phys.addSource("Api.txt"      , "Api");
//     phys.addSource("fractions.txt", "facc");
//
//     if (phys.compute()) {
//         std::cout << phys.value() << " ± " << phys.error() << std::endl;
//     }
// -----------------------------------------------------------------------------
#ifndef PHYSICS_ASYM_MULTI_H
#define PHYSICS_ASYM_MULTI_H

#include <string>
#include <unordered_map>
#include <vector>

struct CorrEntry { double v=0, e=0; };

class PhysicsAsymmetryMulti {
public:
    // expose read-only view of all loaded corrections
    const std::unordered_map<std::string,CorrEntry>&
    corrections() const { return corr_; }

    /// register an input file containing TNamed(<key>,"val,err") or TXT
    void addSource(const std::string& file, const std::string& key);

    /// compute final physics asymmetry using expression passed by user
    ///   expr should be a lambda (const std::unordered_map<std::string,CorrEntry>&)->std::pair<double,double>
    template<typename Expr>
    bool compute(const Expr& expr);

    double value() const { return A_; }
    double error() const { return dA_; }

private:
    CorrEntry loadPair(const std::string& file, const std::string& key);

    std::vector<std::pair<std::string,std::string>> sources_;   // (file,key)
    std::unordered_map<std::string,CorrEntry>       corr_;      // key → {v,e}

    double A_=0, dA_=0;   // result
};

#endif // PHYSICS_ASYM_MULTI_H

