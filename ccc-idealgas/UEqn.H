{
    Info<< "Min U: " << min(mag(U)).value() << ", Max U: " << max(mag(U)).value() << endl;
    Info<< "Min rho: " << min(rho).value() << ", Max rho: " << max(rho).value() << endl;

    dimensionedScalar sigma
    (
        "sigma",
        dimMass/sqr(dimTime),
        0.072 // Surface tension for R32
    );

    UEqn.relax();
    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK()) * fvc::snGrad(alpha_liquid)
                  - ghf * fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
}
