{
    scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
    scalar maxDeltaT = runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    scalar deltaTFact = min(min(maxCo/0.5, 1.0), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " << runTime.deltaTValue() << endl;
}
