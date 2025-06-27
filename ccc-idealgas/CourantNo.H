{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    // Check if mesh is valid
    if (mesh.V().size() == 0)
    {
        FatalErrorInFunction
            << "Invalid mesh: No cells found."
            << abort(FatalError);
    }

    // Check if phi is valid
    if (phi.internalField().size() == 0)
    {
        FatalErrorInFunction
            << "Flux field (phi) has no internal field data."
            << abort(FatalError);
    }

    // Read maxCo from controlDict
    scalar maxCo = runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.1);

    // Compute flux magnitude
    surfaceScalarField sumPhi = mag(phi);

    // Check for zero or negative cell volumes
    scalarField cellVolumes = mesh.V().field();
    if (min(cellVolumes) <= 0.0)
    {
        FatalErrorInFunction
            << "Invalid mesh: Zero or negative cell volumes detected."
            << abort(FatalError);
    }

    // Maximum Courant number with safeguard
    scalarField phiOverV = sumPhi / (cellVolumes + SMALL);
    CoNum = 0.5 * gMax(phiOverV) * runTime.deltaTValue();

    // Mean Courant number with safeguard
    meanCoNum = 0.5 * (gSum(sumPhi) / (gSum(cellVolumes) + SMALL)) * runTime.deltaTValue();

    // Cap Courant number by adjusting time step
    if (CoNum > maxCo && runTime.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        scalar newDeltaT = runTime.deltaTValue() * maxCo / (CoNum + SMALL);
        runTime.setDeltaT(newDeltaT);
        CoNum = maxCo; // Update CoNum to reflect capped value
        Info<< "Adjusted deltaT to " << newDeltaT << " to cap CoNum at " << maxCo << endl;
    }

    // Output
    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;
}
