nCorr = stressControl.getOrDefault<int>("nCorrectors", 1);
convergenceTolerance = stressControl.get<scalar>("D");
