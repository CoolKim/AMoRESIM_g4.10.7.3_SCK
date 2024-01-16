#ifndef __AmoreSimLibGitRevision_H_
#define __AmoreSimLibGitRevision_H_

#include <iostream>

class AmoreSimGitRevision {
  public:
    static void PrintGitInfo();
    static const char *GetRevisionHash();
    static const char *GetBranchName();

  private:
    AmoreSimGitRevision(){}; // Don't create an object of this cless!
    ~AmoreSimGitRevision(){};
    static const char *fgceRevisionHash;
    static const char *fgceBranchName;
    static const char *fgceLibraryName;
};

#endif
