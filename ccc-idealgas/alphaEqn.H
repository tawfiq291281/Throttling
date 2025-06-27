#include "fvc.H"
#include "fvScalarMatrix.H"
#include "MULES.H"

{
    volScalarField& alpha1 = alpha_liquid;
    volScalarField alpha2 = 1.0 - alpha1;

    word alphaScheme("div(phi,alpha.liquid)");
    word alpharScheme("div(phirb,alpha.liquid)");

    Info<< "Min phi: " << min(phi).value() << ", Max phi: " << max(phi).value() << endl;

    // Declare phiAlpha for MULES
    surfaceScalarField phiAlpha
    (
        IOobject
        (
            "phiAlpha",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi
    );

    // Create scalarFields for psiMax and psiMin
    scalarField psiMax(mesh.nCells(), 1.0);
    scalarField psiMin(mesh.nCells(), 0.0);

    // Explicitly limit alpha1 for MULES
    MULES::explicitSolve(rho, alpha1, phi, phiAlpha, psiMax, psiMin);

    volScalarField mDotAlphal = mixture->mDotAlphal();
    Info<< "mDotAlphal dimensions: " << mDotAlphal.dimensions() << endl;
    Info<< "Min mDotAlphal: " << min(mDotAlphal).value() << ", Max mDotAlphal: " << max(mDotAlphal).value() << endl;
    volScalarField Su = mDotAlphal/rho1; // Normalize by liquid density
    volScalarField Sp = -mDotAlphal/rho1; // Normalize by liquid density
    Su = min(max(Su, dimensionedScalar("SuMin", Su.dimensions(), -1e3)), dimensionedScalar("SuMax", Su.dimensions(), 1e3));
    Sp = min(max(Sp, dimensionedScalar("SpMin", Sp.dimensions(), -1e3)), dimensionedScalar("SpMax", Sp.dimensions(), 1e3));
    Info<< "Min Su: " << min(Su).value() << ", Max Su: " << max(Su).value() << endl;
    Info<< "Min Sp: " << min(Sp).value() << ", Max Sp: " << max(Sp).value() << endl;

    for (int aCorr = 0; aCorr < nAlphaCorr; aCorr++)
    {
        surfaceScalarField phic(mag(phi/mesh.magSf()));
        phic = min(interface.cAlpha()*phic, max(phic));
        surfaceScalarField phir(phic*interface.nHatf());

        surfaceScalarField phiAlphaCorr
        (
            fvc::flux(phir, alpha1, alpharScheme)
        );

        phiAlpha += phiAlphaCorr;

        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha1) + fvm::div(phi, alpha1, alphaScheme) + fvc::div(phiAlphaCorr) == Su + fvm::Sp(Sp, alpha1)
        );

        Info<< "alphaEqn matrix diagonal: " << alphaEqn.diagonal() << endl;
        Info<< "alphaEqn matrix symmetric: " << alphaEqn.symmetric() << endl;
        Info<< "alphaEqn matrix asymmetric: " << alphaEqn.asymmetric() << endl;

        dictionary solverDict;
        solverDict.add("solver", "PBiCGStab");
        solverDict.add("preconditioner", "DILU");
        solverDict.add("tolerance", 1e-8); // Tightened
        solverDict.add("relTol", 0.01);
        solverDict.add("maxIter", 2000);

        Info<< "Solving alphaEqn with PBiCGStab, iteration " << aCorr + 1 << endl;
        alphaEqn.solve(solverDict);

        // Explicitly bound alpha1
        alpha1 = min(max(alpha1, 0.0), 1.0);
        alpha2 = 1.0 - alpha1;

        Info<< "After iteration " << aCorr + 1 << ": Min(alpha1) = " << min(alpha1).value() << ", Max(alpha1) = " << max(alpha1).value() << endl;
    }

    // Update density using mixture's rho2
    rho = alpha1*mixture->rho1() + alpha2*mixture->rho2();
    rho = max(rho, dimensionedScalar("rhoMin", rho.dimensions(), 1.0));
    rho = min(rho, dimensionedScalar("rhoMax", rho.dimensions(), 960.0));

    rhoPhi = fvc::interpolate(alpha1*(mixture->rho1() - mixture->rho2()))*phi 
           + fvc::interpolate(mixture->rho2())*phi;

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << " Min(alpha1) = " << min(alpha1).value()
        << " Max(alpha1) = " << max(alpha1).value()
        << " Min(alpha2) = " << min(alpha2).value()
        << " Max(alpha2) = " << max(alpha2).value()
        << endl;
}
