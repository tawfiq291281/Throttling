{
    const volScalarField& alpha1 = mixture->alpha1();
    volScalarField alpha2 = 1.0 - alpha1;

    volScalarField cp = mixture->cp();
    cp = max(cp, dimensionedScalar("cpMin", cp.dimensions(), 500.0));
    cp = min(cp, dimensionedScalar("cpMax", cp.dimensions(), 5000.0));

    dimensionedScalar kappa1
    (
        transportProperties.subDict("liquid").get<dimensionedScalar>("kappa")
    );
    dimensionedScalar kappa2
    (
        transportProperties.subDict("vapour").get<dimensionedScalar>("kappa")
    );
    volScalarField kappa = alpha1*kappa1 + alpha2*kappa2;
    kappa = max(kappa, dimensionedScalar("kappaMin", kappa.dimensions(), 0.01));

    const Lee& lee = dynamic_cast<const Lee&>(*mixture);
    volScalarField heSource = mixture->heSource();
    volScalarField mDotEvaporation = lee.mDotEvaporation()(); // Use accessor
    volScalarField mDotCondensation = lee.mDotCondensation()(); // Use accessor
    volScalarField mDotAlphal = mixture->mDotAlphal(); // For consistency check

    volScalarField h_old = h;

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(rhoPhi, h)
      - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), h)
      - fvm::laplacian(kappa/cp, h)
     ==
        heSource
      + fvOptions(rho, h)
    );

    EEqn.relax();
    fvOptions.constrain(EEqn);
    EEqn.solve();
    fvOptions.correct(h);

    Info<< "Before capping: Min h: " << min(h).value() << ", Max h: " << max(h).value() << endl;
    h = max(h, dimensionedScalar("hMin", h.dimensions(), 200000));
    h = min(h, dimensionedScalar("hMax", h.dimensions(), 400000)); // Tighten for 354,105 J/kg
    h.correctBoundaryConditions();

    T = mixture->T(h);
    Info<< "Before capping: Min T: " << min(T).value() << ", Max T: " << max(T).value() << endl;
    T = max(T, dimensionedScalar("Tmin", T.dimensions(), 272.95)); // Align with 0.908 MPa
    T = min(T, dimensionedScalar("Tmax", T.dimensions(), 320.0));
    Info<< "After capping: Min T: " << min(T).value() << ", Max T: " << max(T).value() << endl;
    T.correctBoundaryConditions();
    if (max(T).value() > 320.0 || min(T).value() < 272.95)
    {
        WarningInFunction << "Temperature out of bounds after capping: Min T: " << min(T).value()
                          << ", Max T: " << max(T).value() << endl;
    }

    const volScalarField& p_rgh = mesh.lookupObject<volScalarField>("p_rgh");
    const dimensionedScalar P0("P0", dimPressure, 1e5);
    volScalarField p_rgh_limited = max(
        min(p_rgh, dimensionedScalar("pMax", p_rgh.dimensions(), 3.2e6)),
        dimensionedScalar("pMin", p_rgh.dimensions(), 9e5)
    );
    volScalarField Tsat = lee.getB() / (lee.getA() - log10(p_rgh_limited / P0)) - lee.getC();
    Info<< "Before capping: Min Tsat: " << min(Tsat).value() << ", Max Tsat: " << max(Tsat).value() << endl;
    Tsat = max(min(Tsat, dimensionedScalar("TsatMax", Tsat.dimensions(), 320.0)), dimensionedScalar("TsatMin", T.dimensions(), 272.95));
    Info<< "After capping: Min Tsat: " << min(Tsat).value() << ", Max Tsat: " << max(Tsat).value() << endl;
    Info<< "Min mDotEvaporation: " << min(mDotEvaporation).value() << ", Max mDotEvaporation: " << max(mDotEvaporation).value() << endl;
    Info<< "Min mDotCondensation: " << min(mDotCondensation).value() << ", Max mDotCondensation: " << max(mDotCondensation).value() << endl;
    Info<< "Min mDotAlphal: " << min(mDotAlphal).value() << ", Max mDotAlphal: " << max(mDotAlphal).value() << endl;
    Info<< "Min heSource: " << min(heSource).value() << ", Max heSource: " << max(heSource).value() << endl;
    Info<< "After capping: Min h: " << min(h).value() << ", Max h: " << max(h).value() << endl;
}
