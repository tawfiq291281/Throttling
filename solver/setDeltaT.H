{
    if (adjustTimeStep)
    {
        scalar maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", 1e-3);
        scalar minDeltaT = runTime.controlDict().lookupOrDefault<scalar>("minDeltaT", 1e-6);
        scalar targetCo = 0.1;

        reduce(maxCo, maxOp<scalar>());

        scalar deltaT = maxDeltaT;
        if (maxCo > targetCo)
        {
            deltaT = min(maxDeltaT, maxDeltaT/(maxCo/targetCo + 0.1));
        }

        deltaT = max(minDeltaT, deltaT);
        runTime.setDeltaT(deltaT);
    }
}
