// File: plot_runbyrun_asym.C
//
// Usage:
//   root -l
//   .x plot_runbyrun_asym.C
//
// or
//   .L plot_runbyrun_asym.C
//   plot_runbyrun_asym("raw_asymmetry_GEN2_He3.txt",
//                      "raw_proton_asymmetry_GEN2_He3.txt");

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>
#include <TStyle.h>
#include <TString.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

#include <map>
#include <set>
#include <TLatex.h>

struct AsymPoint {
    double run  = 0.0;
    double asym = 0.0;   // stored in %
    double err  = 0.0;   // stored in %
};

std::vector<AsymPoint> ReadNeutronFile(const char* filename) {
    std::vector<AsymPoint> pts;

    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Error: could not open neutron file: " << filename << std::endl;
        return pts;
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);

        int run;
        double Np, Nm, Araw, dAraw, beamPol, dBeam, targetPol, dTgt;
        ss >> run >> Np >> Nm >> Araw >> dAraw >> beamPol >> dBeam >> targetPol >> dTgt;

        if (ss.fail()) continue;

        AsymPoint p;
        p.run  = run;
        p.asym = 100.0 * Araw;
        p.err  = 100.0 * dAraw;
        pts.push_back(p);
    }

    return pts;
}

std::vector<AsymPoint> ReadProtonFile(const char* filename) {
    std::vector<AsymPoint> pts;

    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Error: could not open proton file: " << filename << std::endl;
        return pts;
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);

        int run;
        double Np, Nm, Araw, dAraw;
        ss >> run >> Np >> Nm >> Araw >> dAraw;

        if (ss.fail()) continue;

        AsymPoint p;
        p.run  = run;
        p.asym = 100.0 * Araw;
        p.err  = 100.0 * dAraw;
        pts.push_back(p);
    }

    return pts;
}

// TGraphErrors* MakeGraph(const std::vector<AsymPoint>& pts, Color_t color, Style_t markerStyle) {
//     if (pts.empty()) return nullptr;

//     auto* g = new TGraphErrors((int)pts.size());
//     for (int i = 0; i < (int)pts.size(); i++) {
//         g->SetPoint(i, pts[i].run, pts[i].asym);
//         g->SetPointError(i, 0.0, pts[i].err);
//     }

//     g->SetMarkerStyle(markerStyle);
//     g->SetMarkerSize(0.7);
//     g->SetMarkerColor(color);
//     g->SetLineColor(color);
//     g->SetLineWidth(1);

//     return g;
//}

TGraphErrors* MakeGraph(const std::vector<AsymPoint>& pts,
                        const std::map<int,double>& runToPseudo,
                        Color_t color, Style_t markerStyle) {
    if (pts.empty()) return nullptr;

    auto* g = new TGraphErrors((int)pts.size());
    for (int i = 0; i < (int)pts.size(); i++) {
        int run = (int)std::lround(pts[i].run);
        auto it = runToPseudo.find(run);
        if (it == runToPseudo.end()) continue;

        g->SetPoint(i, it->second, pts[i].asym);
        g->SetPointError(i, 0.0, pts[i].err);
    }

    g->SetMarkerStyle(markerStyle);
    g->SetMarkerSize(0.8);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);

    return g;
}

void DrawHorizontalGuides(double xmin, double xmax, double yMin=-10.0, double yMax=10.0, double step=2.0) {
    for (double y = yMin; y <= yMax + 1e-6; y += step) {
        auto* line = new TLine(xmin, y, xmax, y);
        line->SetLineColor(kGray+1);

        if (std::fabs(y) < 1e-9) {
            line->SetLineStyle(1);
            line->SetLineWidth(2);
        } else {
            line->SetLineStyle(2);
            line->SetLineWidth(1);
        }
        line->Draw("same");
    }
}

void plot_runbyrun_asym(const char* neutronFile = "/mnt/data/raw_asymmetry_GEN2_He3.txt",
                        const char* protonFile  = "/mnt/data/raw_proton_asymmetry_GEN2_He3.txt") {

    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(3);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    auto neutron = ReadNeutronFile(neutronFile);
    auto proton  = ReadProtonFile(protonFile);

    if (neutron.empty() && proton.empty()) {
        std::cerr << "No valid points found in either file." << std::endl;
        return;
    }

    // Build sorted unique run list from both files
    std::set<int> uniqueRuns;
    for (const auto& p : neutron) uniqueRuns.insert((int)std::lround(p.run));
    for (const auto& p : proton ) uniqueRuns.insert((int)std::lround(p.run));

    std::vector<int> runList(uniqueRuns.begin(), uniqueRuns.end());

    // Map actual run -> pseudo run index
    std::map<int,double> runToPseudo;
    for (int i = 0; i < (int)runList.size(); i++) {
        runToPseudo[runList[i]] = i + 1;   // start at 1
    }

    double xmin = 0.5;
    double xmax = runList.size() + 0.5;

    auto* gN = MakeGraph(neutron, runToPseudo, kBlue, 20);
    auto* gP = MakeGraph(proton,  runToPseudo, kRed, 20);

    auto* c = new TCanvas("c_runbyrun_asym", "Run-by-run raw asymmetry", 1400, 500);
    c->SetLeftMargin(0.08);
    c->SetRightMargin(0.03);
    c->SetTopMargin(0.05);
    c->SetBottomMargin(0.14);

    auto* frame = c->DrawFrame(xmin, -20.0, xmax, 20.0);
    frame->GetXaxis()->SetTitle("Run number (good runs only)");
    frame->GetYaxis()->SetTitle("Raw Asymmetry (%)");
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();

    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.03);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetNdivisions(510);

    //frame->GetXaxis()->SetTitle("Pseudo run number");
    frame->GetXaxis()->SetLabelSize(0.0);     // hide ROOT's default labels
    frame->GetXaxis()->SetTickLength(0.02);

    int nLabels = 10;
    double yText = -21.5;   // place labels a bit below the frame

    for (int j = 0; j < nLabels; j++) {
        int idx = (int)std::round(j * (runList.size()-1.0) / (nLabels-1.0));
        double x = idx + 1;   // pseudo run position

        TLatex *lab = new TLatex(x, yText, Form("%d", runList[idx]));
        lab->SetTextAlign(23);    // centered
        lab->SetTextSize(0.025);
        lab->Draw();
    }

    DrawHorizontalGuides(xmin, xmax, -10.0, 10.0, 2.0);

    if (gN) gN->Draw("P SAME");
    if (gP) gP->Draw("P SAME");

    auto* leg = new TLegend(0.80, 0.82, 0.97, 0.95);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);
    if (gN) leg->AddEntry(gN, "(e,e'n) events", "p");
    if (gP) leg->AddEntry(gP, "(e,e'p) events", "p");
    leg->Draw();

    c->SaveAs("runbyrun_raw_asymmetry_GEN2_He3.pdf");
    c->SaveAs("runbyrun_raw_asymmetry_GEN2_He3.png");
}
