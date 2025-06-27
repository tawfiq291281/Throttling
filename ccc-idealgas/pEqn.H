{
      volScalarField rAU = 1.0/UEqn.A();
      rAU = max(rAU, dimensionedScalar("rAUMin", rAU.dimensions(), 1e-6));
      rAU = min(rAU, dimensionedScalar("rAUMax", rAU.dimensions(), 1e10)); // Tightened
      Info<< "Min rAU: " << min(rAU).value() << ", Max rAU: " << max(rAU).value() << endl;

      surfaceScalarField rAUf = fvc::interpolate(rAU);

      volVectorField HbyA = rAU*UEqn.H();

      surfaceScalarField phiHbyA
      (
          IOobject
          (
              "phiHbyA",
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          fvc::flux(HbyA)
      );

      adjustPhi(phiHbyA, U, p_rgh);

      // Bound p_rgh before solving
      p_rgh = max(p_rgh, dimensionedScalar("p_rghMin", p_rgh.dimensions(), 9e5)); // Align with outlet
      p_rgh = min(p_rgh, dimensionedScalar("p_rghMax", p_rgh.dimensions(), 3.2e6));
      p_rgh.correctBoundaryConditions();

      while (pimple.correctNonOrthogonal())
      {
          fvScalarMatrix p_rghEqn
          (
              fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
          );

          p_rghEqn.setReference(pRefCell, pRefValue);
          p_rghEqn.solve(mesh.solver("p_rgh"));

          if (pimple.finalNonOrthogonalIter())
          {
              phi = phiHbyA - p_rghEqn.flux();
              phi = max(phi, dimensionedScalar("phiMin", phi.dimensions(), -1e-3)); // Tightened
              phi = min(phi, dimensionedScalar("phiMax", phi.dimensions(), 1e-3));
              rhoPhi = fvc::interpolate(rho)*phi;
              rhoPhi = max(rhoPhi, dimensionedScalar("rhoPhiMin", rhoPhi.dimensions(), -960*1e-3));
              rhoPhi = min(rhoPhi, dimensionedScalar("rhoPhiMax", rhoPhi.dimensions(), 960*1e-3));
          }
      }

      U = HbyA - rAU*fvc::grad(p_rgh);
      U.correctBoundaryConditions();
      fvOptions.correct(U);

      p = max(p, dimensionedScalar("pMin", p.dimensions(), 9e5)); // Align with outlet
      p = min(p, dimensionedScalar("pMax", p.dimensions(), 3.2e6));
      p.correctBoundaryConditions();

      rho2 = mixture->rho2();
      rho2 = max(rho2, dimensionedScalar("rho2Min", rho2.dimensions(), 50.0)); // Align with R32 vapor
      rho2 = min(rho2, dimensionedScalar("rho2Max", rho2.dimensions(), 100.0));
      rho = alpha_liquid*rho1 + alpha_vapour*rho2;
      rho = max(rho, dimensionedScalar("rhoMin", rho.dimensions(), 1.0));
      rho = min(rho, dimensionedScalar("rhoMax", rho.dimensions(), 960.0));
      rho.correctBoundaryConditions();

      Info<< "Min p_rgh: " << min(p_rgh).value() << ", Max p_rgh: " << max(p_rgh).value() << endl;
      Info<< "Min p: " << min(p).value() << ", Max p: " << max(p).value() << endl;
      Info<< "Min rho: " << min(rho).value() << ", Max rho: " << max(rho).value() << endl;
}
