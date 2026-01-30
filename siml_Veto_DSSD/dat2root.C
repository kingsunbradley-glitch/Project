#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "TFile.h"
#include "TTree.h"

void dat2root() {
   
    std::string inputFileName = "siml_Veto_DSSD.dat";
    std::string outputFileName = "Veto_DSSD.root";
    std::string treeName = "siml";

    std::ifstream infile(inputFileName.c_str());
    if (!infile.is_open()) {
        std::cout << "Error: Cannot open file " << inputFileName << std::endl;
        return;
    }

   
    std::string headerLine;
    if (!std::getline(infile, headerLine)) {
        std::cout << "Error: File is empty." << std::endl;
        return;
    }

    std::stringstream ss(headerLine);
    std::string colName;
    std::vector<std::string> headers;
    while (ss >> colName) {
        headers.push_back(colName);
    }
    
    int numColumns = headers.size();
    std::cout << "Found " << numColumns << " columns in header." << std::endl;

    
    TFile *outfile = new TFile(outputFileName.c_str(), "RECREATE");
    TTree *tree = new TTree(treeName.c_str(), "Data from siml_Veto_DSSD.dat");

   
    Int_t idx;
    std::vector<Double_t> values(numColumns - 1);
    tree->Branch(headers[0].c_str(), &idx, (headers[0] + "/I").c_str());

    
    for (int i = 1; i < numColumns; ++i) {
        std::string leafType = headers[i] + "/D";
        tree->Branch(headers[i].c_str(), &values[i-1], leafType.c_str());
    }

    //读取数据
    while (infile >> idx) {
        for (int i = 0; i < numColumns - 1; ++i) {
            infile >> values[i];
        }
        tree->Fill();
    }
    tree->Write();
    outfile->Close();
    infile.close();

    std::cout << "Successfully created " << outputFileName << " with tree '" << treeName << "'" << std::endl;
}