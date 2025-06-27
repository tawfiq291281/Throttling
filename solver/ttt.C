#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "myPhaseChangeTwoPhaseMixture.H"
#include "Lee.H"
#include "CMULES.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "OFstream.H"
#include "cpuTime.H"
#include "PstreamReduceOps.H"
#include "VectorTensorOps.H"
#include "alphaControlsVars.H"

// Critical parameters
label nAlphaCorr = 2;
label nAlphaSubCycles = 4;
bool adjustTimeStep = true;
scalar maxCo = 0.02; // Tightened for stability

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Parallel solver for two fluids with phase-change,\n"
        "constant enthalpy throttling for R-32 in a 3D spiral tube."
    );

    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    if (Pstream::parRun())
        Info<< "Running in parallel on " << Pstream::nProcs() << " processors" << endl;
    else
        Info<< "Running in serial mode" << endl;

    #include "createControl.H"

    IOdictionary fvSchemes
    (
        IOobject
        (
            "fvSchemes",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    {
        dictionary alphaControls = mesh.solutionDict().subOrEmptyDict("alpha");
        nAlphaCorr = alphaControls.getOrDefault<label>("nAlphaCorr", 2);
        nAlphaSubCycles = alphaControls.getOrDefault<label>("nAlphaSubCycles", 4);
    }

    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    if (!turbulence.valid())
        FatalErrorInFunction << "Turbulence model not initialized correctly." << exit(FatalError);

    cpuTime executionTime;
    scalar elapsedCpuTime = 0.0;

    dimensionedScalar maxUmag("maxUmag", dimVelocity, 15.0);

    label pRefCell = 0;
    scalar pRefValue = 908000; // Outlet pressure
    IOdictionary fvSolution
    (
        IOobject
        (
            "fvSolution",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    if (fvSolution.subOrEmptyDict("PIMPLE").found("pRefCell"))
    {
        pRefCell = fvSolution.subDict("PIMPLE").get<label>("pRefCell");
        pRefValue = fvSolution.subDict("PIMPLE").get<scalar>("pRefValue");
    }

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"

        if (adjustTimeStep)
        {
            scalar minDeltaT = 1e-8;
            scalar maxDeltaT = 1e-9; // Tightened
            scalar targetCo = 0.02; // Tightened

            reduce(maxCo, maxOp<scalar>());

            scalar deltaT = maxDeltaT;
            if (maxCo > targetCo)
                deltaT = min(maxDeltaT, maxDeltaT/(maxCo/targetCo + 0.1));

            deltaT = max(minDeltaT, deltaT);
            runTime.setDeltaT(deltaT, false);
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
            mixture->correct();
            #include "alphaEqn.H"
            interface.correct();
            #include "EEqn.H"

            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvc::div(rhoPhi, U)
              + turbulence->divDevRhoReff(rho,U)
             ==
                fvOptions(rho, U)
              - rho*g
              + fvc::reconstruct(fvc::interpolate(interface.sigmaK()) * fvc::snGrad(alpha_liquid) * mesh.magSf())
            );

            #include "UEqn.H"
            #include "pEqn.H"

            if (pimple.turbCorr())
            {
                Info<< "Executing turbulence correction\n" << endl;
                turbulence->correct();
                turbProd = max(turbulence->nuEff() * magSqr(symm(fvc::grad(U))), 
                               dimensionedScalar("small", dimLength*dimLength/(dimTime*dimTime*dimTime), 1e-10));
            }

            // Cap velocity magnitude
            volScalarField UMag = mag(U);
            dimensionedScalar UMagMax("UMagMax", dimVelocity, 100);
            volScalarField scale = min(UMagMax / max(UMag, dimensionedScalar("smallU", dimVelocity, 1e-6)), 1.0);
            U *= scale;
            U.correctBoundaryConditions();
            fvOptions.correct(U);

            Info<< "After PIMPLE loop: Min U: " << min(mag(U)).value() << ", Max U: " << max(mag(U)).value() << endl;
        }

        elapsedCpuTime = executionTime.elapsedCpuTime();
        reduce(elapsedCpuTime, sumOp<scalar>());
        elapsedCpuTime /= Pstream::nProcs();

        if (Pstream::master())
            Info<< "ExecutionTime = " << elapsedCpuTime << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl;

        runTime.write();

        if (Pstream::master())
            runTime.printExecutionTime(Info);
    }

    if (Pstream::master())
        Info<< "End\n" << nl << endl;

    return 0;
}
