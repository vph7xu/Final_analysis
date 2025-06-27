#include "BranchVars.h"
#include "CutConfig.h"
#include "CutManager.h"
#include "RawAsymmetry.h"

#include <TChain.h>
#include <iostream>

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <root-file(s)> <cut-file.json>\n";
        return 1;
    }

    /* 1. build the chain */
    TChain ch("Tout");                    // change tree name if needed
    for (int i = 1; i < argc - 1; ++i) ch.Add(argv[i]);

    /* 2. attach branches */
    BranchVars v;
    v.attach(&ch);

    /* 3. cuts */
    CutConfig cfg;
    if (!cfg.load(argv[argc - 1])) {
        std::cerr << "Cannot read cuts file " << argv[argc - 1] << '\n';
        return 2;
    }
    CutManager cut(cfg);

    /* 4. run the module */
    RawAsymmetry mod(cut);
    mod.process(ch, v);
}

