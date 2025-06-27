#include "fvCFD.H"
#include "myPhaseChangeTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "interfaceProperties.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field h\n" << endl;
volScalarField h
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating face flux field phi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::flux(U)
);

Info<< "Reading field alpha.liquid\n" << endl;
volScalarField alpha_liquid
(
    IOobject
    (
        "alpha.liquid",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField alpha_vapour
(
    IOobject
    (
        "alpha.vapour",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1.0 - alpha_liquid
);

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Reading specific gas constant for vapour\n" << endl;
dimensionedScalar R
(
    transportProperties.subDict("vapour").get<dimensionedScalar>("R")
);

Info<< "Creating phase change mixture\n" << endl;
autoPtr<myPhaseChangeTwoPhaseMixture> mixture = myPhaseChangeTwoPhaseMixture::New(U, phi, alpha_liquid, alpha_vapour);

Info<< "Creating interface properties\n" << endl;
interfaceProperties interface(alpha_liquid, U, mixture());

// Liquid density (constant for R-32 liquid)
dimensionedScalar rho1
(
    transportProperties.subDict("liquid").get<dimensionedScalar>("rho")
);

// Check rho1 dimensions
if (rho1.dimensions() != dimensionSet(1, -3, 0, 0, 0, 0, 0))
{
    FatalErrorInFunction
        << "Invalid dimensions for rho1: " << rho1.dimensions()
        << ". Expected [1 -3 0 0 0 0 0]." << exit(FatalError);
}

// Declare p before rho2
Info<< "Calculating field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh // Initialize p with p_rgh
);

// Vapor density from mixture
volScalarField rho2 = mixture->rho2();

// Check rho2 dimensions
if (rho2.dimensions() != dimensionSet(1, -3, 0, 0, 0, 0, 0))
{
    FatalErrorInFunction
        << "Invalid dimensions for rho2: " << rho2.dimensions()
        << ". Expected [1 -3 0 0 0 0 0]." << exit(FatalError);
}

// Total density
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha_liquid*rho1 + alpha_vapour*rho2
);

// Check rho dimensions
if (rho.dimensions() != dimensionSet(1, -3, 0, 0, 0, 0, 0))
{
    FatalErrorInFunction
        << "Invalid dimensions for rho: " << rho.dimensions()
        << ". Expected [1 -3 0 0 0 0 0]." << exit(FatalError);
}

Info<< "Calculating field nu\n" << endl;
volScalarField nu
(
    IOobject
    (
        "nu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mixture->nu()
);

Info<< "Calculating field rhoPhi\n" << endl;
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(rho) * phi
);

Info<< "Reading gravitational acceleration\n" << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Calculating field gh\n" << endl;
volScalarField gh("gh", g & mesh.C());

Info<< "Calculating field ghf\n" << endl;
surfaceScalarField ghf("ghf", fvc::interpolate(gh));

// Update p
p = p_rgh + rho*gh;

Info<< "Creating turbulence model (k-omega)\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New
    (
        U,
        phi,
        mixture()
    )
);

Info<< "Reading field k\n" << endl;
volScalarField k
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field omega\n" << endl;
volScalarField omega
(
    IOobject
    (
        "omega",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating field nut\n" << endl;
volScalarField nut
(
    IOobject
    (
        "nut",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->nut()
);

Info<< "Creating fvOptions\n" << endl;
fv::options fvOptions(mesh);

Info<< "Creating turbulent production field\n" << endl;
volScalarField turbProd
(
    IOobject
    (
        "turbProd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    max(turbulence->nuEff() * magSqr(symm(fvc::grad(U))), dimensionedScalar("small", dimensionSet(0, 2, -3, 0, 0, 0, 0), 1e-6))
);
