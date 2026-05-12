#include "createControl.H"

const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");

int nCorr = stressControl.getOrDefault<int>("nCorrectors", 1);

scalar convergenceTolerance(stressControl.get<scalar>("D"));
