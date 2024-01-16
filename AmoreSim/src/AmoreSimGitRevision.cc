#include <unistd.h>

#include "AmoreSim/AmoreSimGitRevision.hh"

#define STRINGFY(X) #X
#define TOSTRING(X) STRINGFY(X)

const char *AmoreSimGitRevision::fgceRevisionHash = TOSTRING(AmoreSim_GIT_COMMIT_HASH);
const char *AmoreSimGitRevision::fgceBranchName   = TOSTRING(AmoreSim_GIT_BRANCH);
const char *AmoreSimGitRevision::fgceLibraryName  = "AmoreSim library";

const char *AmoreSimGitRevision::GetRevisionHash() { return fgceRevisionHash; }
const char *AmoreSimGitRevision::GetBranchName() { return fgceBranchName; }
void AmoreSimGitRevision::PrintGitInfo() {
    std::cout << fgceLibraryName << ": <Branch: " << fgceBranchName
              << ", Revision hash: " << fgceRevisionHash << ">" << std::endl;
}

int main() {
    using namespace std;
    cout << "       AmoreSim library for MC simulation of AMoRE project" << endl;
    cout << "       Author: Center for Underground Physics (CUP), Korea" << endl;
    cout << "=================================================================" << endl;
    cout << "This library provides classes of Geant4 for a simulation of AMoRE" << endl
         << "(Advenced 100Mo-based Rare process Experiment)." << endl
         << "Classes in this library is an implementation of functionality for" << endl
         << "AMoRE-specific geometries, routines for information recording and" << endl
         << "Geant4 messenger commands." << endl;
    cout << "=================================================================" << endl;
    cout << "AmoreSim library: "
         << "<Branch: " << TOSTRING(AmoreSim_GIT_BRANCH)
         << ", Revision hash: " << TOSTRING(AmoreSim_GIT_COMMIT_HASH) << ">" << endl;
    _exit(0);
}
