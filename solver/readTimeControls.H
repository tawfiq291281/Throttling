{
    adjustTimeStep = runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

    // scalar maxDeltaT =
    //     runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
}
