#include "Lee.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

namespace Foam
{
    defineTypeNameAndDebug(Lee, 0);
    addToRunTimeSelectionTable
    (
        myPhaseChangeTwoPhaseMixture,
        Lee,
        components
    );

    Lee::Lee
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& alpha1,
        const volScalarField& alpha2
    )
    :
        myPhaseChangeTwoPhaseMixture(U, phi, alpha1, alpha2),
        r_(dimensionedScalar("r", dimensionSet(0, 0, -1, 0, 0, 0, 0), 2.0)), // Increased for stronger evaporation
        L_(dimensionedScalar("L", dimensionSet(0, 2, -2, 0, 0, 0, 0), 219650.0)), // R32 latent heat at 0.908 MPa
        a_(phaseChangeDict_.get<dimensionedScalar>("a")),
        b_(phaseChangeDict_.get<dimensionedScalar>("b")),
        c_(phaseChangeDict_.get<dimensionedScalar>("c")),
        T_
        (
            IOobject
            (
                "T",
                U_.time().timeName(),
                U_.mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            U_.mesh()
        )
    {
        // Verify dimensions of input parameters
        if (r_.dimensions() != dimensionSet(0, 0, -1, 0, 0, 0, 0))
            FatalErrorInFunction << "r has wrong dimensions: " << r_.dimensions() 
                                << ". Expected [0 0 -1 0 0 0 0]." << exit(FatalError);
        
        if (L_.dimensions() != dimensionSet(0, 2, -2, 0, 0, 0, 0))
            FatalErrorInFunction << "L has wrong dimensions: " << L_.dimensions() 
                                << ". Expected [0 2 -2 0 0 0 0]." << exit(FatalError);
        
        if (a_.dimensions() != dimless)
            FatalErrorInFunction << "a has wrong dimensions: " << a_.dimensions() 
                                << ". Expected [0 0 0 0 0 0 0]." << exit(FatalError);
        
        if (b_.dimensions() != dimTemperature)
            FatalErrorInFunction << "b has wrong dimensions: " << b_.dimensions() 
                                << ". Expected [0 0 0 1 0 0 0]." << exit(FatalError);
        
        if (c_.dimensions() != dimTemperature)
            FatalErrorInFunction << "c has wrong dimensions: " << c_.dimensions() 
                                << ". Expected [0 0 0 1 0 0 0]." << exit(FatalError);
        
        Info<< "Lee::Lee: r_ dimensions: " << r_.dimensions() << ", value: " << r_.value() << endl;
        Info<< "Lee::Lee: L_ dimensions: " << L_.dimensions() << ", value: " << L_.value() << endl;
        Info<< "Lee::Lee: a_ dimensions: " << a_.dimensions() << ", value: " << a_.value() << endl;
        Info<< "Lee::Lee: b_ dimensions: " << b_.dimensions() << ", value: " << b_.value() << endl;
        Info<< "Lee::Lee: c_ dimensions: " << c_.dimensions() << ", value: " << c_.value() << endl;
    }

    const volScalarField& Lee::T() const
    {
        return T_;
    }

    tmp<volScalarField> Lee::heSource() const
    {
        tmp<volScalarField> mDot = mDotAlphal();
        tmp<volScalarField> tHeSource = L_ * mDot;
        tHeSource.ref() = min(max(tHeSource(), dimensionedScalar("heSourceMin", tHeSource().dimensions(), -1e8)), 
                              dimensionedScalar("heSourceMax", tHeSource().dimensions(), 1e8)); // Increased cap
        Info<< "Min heSource: " << min(tHeSource()).value() 
            << ", Max heSource: " << max(tHeSource()).value() << endl;
        return tHeSource;
    }

    tmp<volScalarField> Lee::T(const volScalarField& h) const
    {
        volScalarField cp = this->cp();
        volScalarField cpSafe = max(cp, dimensionedScalar("small", cp.dimensions(), 1e-6));
        
        if (L_.dimensions() != dimensionSet(0, 2, -2, 0, 0, 0, 0))
            FatalErrorInFunction << "L has wrong dimensions: " << L_.dimensions() 
                                << ". Expected [0 2 -2 0 0 0 0]." << exit(FatalError);
        
        tmp<volScalarField> tT = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "T",
                    U_.time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                (h - L_*alpha2()) / cpSafe
            )
        );
        tT.ref() = max(tT(), dimensionedScalar("Tmin", dimTemperature, 272.95)); // Align with 0.908 MPa
        tT.ref() = min(tT(), dimensionedScalar("Tmax", dimTemperature, 320.0));
        return tT;
    }

    tmp<volScalarField> Lee::mDotAlphal() const
    {
        const volScalarField& h = this->db().lookupObject<volScalarField>("h");
        tmp<volScalarField> tT = this->T(h);
        const volScalarField& T = tT();
        const volScalarField& p_rgh = this->db().lookupObject<volScalarField>("p_rgh");
        
        volScalarField p_rgh_limited = max(
            min(p_rgh, dimensionedScalar("pMax", p_rgh.dimensions(), 3.2e6)), // Inlet pressure
            dimensionedScalar("pMin", p_rgh.dimensions(), 9e5) // Outlet pressure
        );
        
        const dimensionedScalar P0("P0", dimPressure, 1e5);
        volScalarField logTerm = log10(p_rgh_limited / P0);
        Info<< "Min log10(p_rgh_limited/P0): " << min(logTerm).value() 
            << ", Max log10(p_rgh_limited/P0): " << max(logTerm).value() << endl;
        
        volScalarField Tsat = b_ / (a_ - logTerm) - c_;
        
        if (Tsat.dimensions() != dimTemperature)
            FatalErrorInFunction << "Tsat has wrong dimensions: " << Tsat.dimensions() 
                                << ". Expected [0 0 0 1 0 0 0]." << exit(FatalError);
        
        Info<< "Before capping: Min Tsat: " << min(Tsat).value() 
            << ", Max Tsat: " << max(Tsat).value() << endl;
        
        dimensionedScalar TsatMin("TsatMin", dimTemperature, 272.95); // Align with 0.908 MPa
        dimensionedScalar TsatMax("TsatMax", dimTemperature, 320.0);
        Tsat = max(min(Tsat, TsatMax), TsatMin);
        
        Info<< "After capping: Min Tsat: " << min(Tsat).value() 
            << ", Max Tsat: " << max(Tsat).value() << endl;
        
        volScalarField deltaT = T - Tsat;
        
        Info<< "Min T-Tsat: " << min(deltaT).value() 
            << ", Max T-Tsat: " << max(deltaT).value() << endl;
        
        volScalarField mDotEvaporation
        (
            IOobject
            (
                "mDotEvaporation",
                U_.time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar("mDotEvaporation", dimDensity/dimTime, 0.0)
        );
        
        volScalarField mDotCondensation
        (
            IOobject
            (
                "mDotCondensation",
                U_.time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar("mDotCondensation", dimDensity/dimTime, 0.0)
        );

        mDotEvaporation = r_ * alpha1() * this->rho1() * pos(T - Tsat) * (T - Tsat) / Tsat;
        mDotCondensation = r_ * alpha2() * this->rho2() * pos(Tsat - T) * (Tsat - T) / Tsat;

        if (mDotEvaporation.dimensions() != dimensionSet(1, -3, -1, 0, 0, 0, 0))
            FatalErrorInFunction << "mDotEvaporation has wrong dimensions: " << mDotEvaporation.dimensions() 
                                << ". Expected [1 -3 -1 0 0 0 0]." << exit(FatalError);
        
        if (mDotCondensation.dimensions() != dimensionSet(1, -3, -1, 0, 0, 0, 0))
            FatalErrorInFunction << "mDotCondensation has wrong dimensions: " << mDotCondensation.dimensions() 
                                << ". Expected [1 -3 -1 0 0 0 0]." << exit(FatalError);

        tmp<volScalarField> mDot = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "mDotAlphal",
                    U_.time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mDotEvaporation - mDotCondensation
            )
        );
        
        mDot.ref() = min(max(mDot(), dimensionedScalar("mDotMin", mDot().dimensions(), -1e3)), 
                         dimensionedScalar("mDotMax", mDot().dimensions(), 1e3)); // Increased cap
        Info<< "Min mDotEvaporation: " << min(mDotEvaporation).value() 
            << ", Max mDotEvaporation: " << max(mDotEvaporation).value() << endl;
        Info<< "Min mDotCondensation: " << min(mDotCondensation).value() 
            << ", Max mDotCondensation: " << max(mDotCondensation).value() << endl;
        Info<< "Min mDotAlphal: " << min(mDot()).value() 
            << ", Max mDotAlphal: " << max(mDot()).value() << endl;
        
        return mDot;
    }

    tmp<volScalarField> Lee::mDotEvaporation() const
    {
        const volScalarField& h = this->db().lookupObject<volScalarField>("h");
        tmp<volScalarField> tT = this->T(h);
        const volScalarField& T = tT();
        const volScalarField& p_rgh = this->db().lookupObject<volScalarField>("p_rgh");
        
        volScalarField p_rgh_limited = max(
            min(p_rgh, dimensionedScalar("pMax", p_rgh.dimensions(), 3.2e6)),
            dimensionedScalar("pMin", p_rgh.dimensions(), 9e5)
        );
        
        const dimensionedScalar P0("P0", dimPressure, 1e5);
        volScalarField Tsat = b_ / (a_ - log10(p_rgh_limited / P0)) - c_;
        Tsat = max(min(Tsat, dimensionedScalar("TsatMax", dimTemperature, 320.0)), 
                   dimensionedScalar("TsatMin", dimTemperature, 272.95));
        
        tmp<volScalarField> tMDotEvap = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "mDotEvaporation",
                    U_.time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                r_ * alpha1() * this->rho1() * pos(T - Tsat) * (T - Tsat) / Tsat
            )
        );
        
        tMDotEvap.ref() = min(max(tMDotEvap(), dimensionedScalar("mDotMin", tMDotEvap().dimensions(), 0.0)), 
                              dimensionedScalar("mDotMax", tMDotEvap().dimensions(), 1e3));
        return tMDotEvap;
    }

    tmp<volScalarField> Lee::mDotCondensation() const
    {
        const volScalarField& h = this->db().lookupObject<volScalarField>("h");
        tmp<volScalarField> tT = this->T(h);
        const volScalarField& T = tT();
        const volScalarField& p_rgh = this->db().lookupObject<volScalarField>("p_rgh");
        
        volScalarField p_rgh_limited = max(
            min(p_rgh, dimensionedScalar("pMax", p_rgh.dimensions(), 3.2e6)),
            dimensionedScalar("pMin", p_rgh.dimensions(), 9e5)
        );
        
        const dimensionedScalar P0("P0", dimPressure, 1e5);
        volScalarField Tsat = b_ / (a_ - log10(p_rgh_limited / P0)) - c_;
        Tsat = max(min(Tsat, dimensionedScalar("TsatMax", dimTemperature, 320.0)), 
                   dimensionedScalar("TsatMin", dimTemperature, 272.95));
        
        tmp<volScalarField> tMDotCond = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "mDotCondensation",
                    U_.time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                r_ * alpha2() * this->rho2() * pos(Tsat - T) * (Tsat - T) / Tsat
            )
        );
        
        tMDotCond.ref() = min(max(tMDotCond(), dimensionedScalar("mDotMin", tMDotCond().dimensions(), 0.0)), 
                              dimensionedScalar("mDotMax", tMDotCond().dimensions(), 1e3));
        return tMDotCond;
    }

    tmp<volScalarField> Lee::vDotP() const
    {
        volScalarField zeroField
        (
            IOobject
            (
                "vDotP",
                U_.time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar(dimDensity/dimTime, Zero)
        );
        return tmp<volScalarField>
        (
            new volScalarField(zeroField)
        );
    }

    tmp<volScalarField> Lee::rho2() const
    {
        const dimensionedScalar R("R", dimensionSet(1, 2, -2, -1, -1, 0, 0), 8.314462618);
        const dimensionedScalar M("M", dimMass/dimMoles, 0.05202); // R32 molar mass
        const dimensionedScalar Tc("Tc", dimTemperature, 351.255); // R32 critical temperature
        const dimensionedScalar Pc("Pc", dimPressure, 5.782e6); // R32 critical pressure
        const scalar omega = 0.276; // R32 acentric factor
        const scalar kappa = 0.37464 + 1.54226*omega - 0.26992*omega*omega;

        const dimensionedScalar a_PR
        (
            "a_PR",
            dimensionSet(1, 5, -2, 0, -2, 0, 0),
            0.45724 * (sqr(R * Tc) / Pc).value()
        );
        
        const dimensionedScalar b_PR
        (
            "b_PR",
            dimensionSet(0, 3, 0, 0, -1, 0, 0),
            0.07780 * ((R * Tc) / Pc).value()
        );

        const volScalarField& p = this->db().lookupObject<volScalarField>("p");
        const volScalarField& T = this->T();
        
        volScalarField p_limited = max(p, dimensionedScalar("pMin", dimPressure, 9e5)); // Outlet pressure
        p_limited = min(p_limited, dimensionedScalar("pMax", dimPressure, 3.2e6)); // Inlet pressure
        
        volScalarField T_safe = max(T, dimensionedScalar("smallT", dimTemperature, 272.95));

        volScalarField Tr = T_safe / Tc;
        volScalarField pr = p_limited / Pc;

        volScalarField alpha = sqr(1.0 + kappa * (1.0 - sqrt(Tr)));

        volScalarField A = a_PR * p_limited / sqr(R * T_safe);
        volScalarField B = b_PR * p_limited / (R * T_safe);

        volScalarField Z
        (
            IOobject
            (
                "Z",
                U_.time().timeName(),
                U_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            U_.mesh(),
            dimensionedScalar("Zinit", dimless, 1.0)
        );

        for (int i = 0; i < 10; ++i)
        {
            volScalarField fZ = Z*Z*Z - (1.0 - B)*Z*Z + (A - 2.0*B - 3.0*sqr(B))*Z - (A*B - sqr(B) - B*sqr(B));
            volScalarField dfZ = 3.0*Z*Z - 2.0*(1.0 - B)*Z + (A - 2.0*B - 3.0*sqr(B));
            Z -= fZ / max(dfZ, dimensionedScalar("small", dimless, 1e-6));
            Z = max(Z, dimensionedScalar("Zinit", dimless, 0.5));
        }

        volScalarField rho_mol = p_limited / (Z * R * T_safe);
        volScalarField rho_mass = rho_mol * M;

        tmp<volScalarField> tRho2 = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho2",
                    U_.time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rho_mass,
                "zeroGradient"
            )
        );
        
        if (tRho2().dimensions() != dimensionSet(1, -3, 0, 0, 0, 0, 0))
        {
            FatalErrorInFunction
                << "Invalid dimensions for tRho2 in Lee::rho2(): " << tRho2().dimensions()
                << ". Expected [1 -3 0 0 0 0 0]." << exit(FatalError);
        }

        Info<< "Lee::rho2() dimensions: " << tRho2().dimensions() << endl;
        
        tRho2.ref() = max(tRho2(), dimensionedScalar("rho2Min", dimDensity, 50.0)); // Align with R32 vapor
        tRho2.ref() = min(tRho2(), dimensionedScalar("rho2Max", dimDensity, 100.0));
        
        Info<< "Min rho2: " << min(tRho2()).value() 
            << ", Max rho2: " << max(tRho2()).value() << endl;
        
        return tRho2;
    }

    void Lee::correct()
    {
        const volScalarField& h = this->db().lookupObject<volScalarField>("h");
        T_ = T(h)();
        T_ = max(T_, dimensionedScalar("Tmin", dimTemperature, 272.95)); // Align with 0.908 MPa
        T_ = min(T_, dimensionedScalar("Tmax", dimTemperature, 320.0));
        myPhaseChangeTwoPhaseMixture::correct();
    }

    bool Lee::read()
    {
        if (myPhaseChangeTwoPhaseMixture::read())
        {
            phaseChangeDict_.readEntry("r", r_);
            phaseChangeDict_.readEntry("a", a_);
            phaseChangeDict_.readEntry("b", b_);
            phaseChangeDict_.readEntry("c", c_);
            
            if (r_.dimensions() != dimensionSet(0, 0, -1, 0, 0, 0, 0))
                FatalErrorInFunction << "r has wrong dimensions: " << r_.dimensions() 
                                    << ". Expected [0 0 -1 0 0 0 0]." << exit(FatalError);
            
            if (a_.dimensions() != dimless)
                FatalErrorInFunction << "a has wrong dimensions: " << a_.dimensions() 
                                    << ". Expected [0 0 0 0 0 0 0]." << exit(FatalError);
            
            if (b_.dimensions() != dimTemperature)
                FatalErrorInFunction << "b has wrong dimensions: " << b_.dimensions() 
                                    << ". Expected [0 0 0 1 0 0 0]." << exit(FatalError);
            
            if (c_.dimensions() != dimTemperature)
                FatalErrorInFunction << "c has wrong dimensions: " << c_.dimensions() 
                                    << ". Expected [0 0 0 1 0 0 0]." << exit(FatalError);
            
            Info<< "Lee::read: r_ dimensions: " << r_.dimensions() << ", value: " << r_.value() << endl;
            Info<< "Lee::read: a_ dimensions: " << a_.dimensions() << ", value: " << a_.value() << endl;
            Info<< "Lee::read: b_ dimensions: " << b_.dimensions() << ", value: " << b_.value() << endl;
            Info<< "Lee::read: c_ dimensions: " << c_.dimensions() << ", value: " << c_.value() << endl;
            return true;
        }
        return false;
    }
} // namespace Foam
