{
    volScalarField contErr(fvc::div(phi));

    dimensionedScalar totalMass = fvc::domainIntegrate(rho);

    scalar sumLocalContErr = fvc::domainIntegrate(mag(contErr)).value()/max(totalMass.value(), 1e-10);
    scalar globalContErr = fvc::domainIntegrate(contErr).value()/max(totalMass.value(), 1e-10);
    cumulativeContErr += globalContErr;

    Info<< "Continuity error: sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr << endl;
}
