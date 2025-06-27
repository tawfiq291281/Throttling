#ifndef myPhaseChangeTwoPhaseMixture_H
#define myPhaseChangeTwoPhaseMixture_H

#include "incompressibleTwoPhaseMixture.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "dimensionedScalar.H"

namespace Foam
{
class myPhaseChangeTwoPhaseMixture
:
    public incompressibleTwoPhaseMixture
{
protected:
    // Critical parameter: Latent heat of phase change
    dimensionedScalar L_;
    // Critical parameters: Liquid and vapor specific heat capacities
    dimensionedScalar cpl1_;
    dimensionedScalar cpv2_;
    // References to phase fraction fields
    const volScalarField& alpha1_;
    const volScalarField& alpha2_;
    // Dictionary for phase change parameters
    dictionary phaseChangeDict_;

public:
    TypeName("myPhaseChangeTwoPhaseMixture");

    declareRunTimeSelectionTable
    (
        autoPtr,
        myPhaseChangeTwoPhaseMixture,
        components,
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& alpha1,
            const volScalarField& alpha2
        ),
        (U, phi, alpha1, alpha2)
    );

    myPhaseChangeTwoPhaseMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& alpha1,
        const volScalarField& alpha2
    );

    static autoPtr<myPhaseChangeTwoPhaseMixture> New
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& alpha1,
        const volScalarField& alpha2
    );

    virtual tmp<volScalarField> heSource() const = 0;
    virtual tmp<volScalarField> T(const volScalarField& h) const = 0;
    virtual tmp<volScalarField> mDotAlphal() const = 0;
    virtual tmp<volScalarField> vDotP() const = 0;
    virtual tmp<volScalarField> rho2() const = 0;

    const dimensionedScalar& L() const { return L_; }
    const dimensionedScalar& cpl1() const { return cpl1_; }
    const dimensionedScalar& cpv2() const { return cpv2_; }

    tmp<volScalarField> cp() const
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "cp",
                    U_.time().timeName(),
                    U_.mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                alpha1_*cpl1_ + alpha2_*cpv2_
            )
        );
    }

    virtual const volScalarField& T() const = 0;

    virtual void correct() {}
    virtual bool read() { return true; }

    virtual ~myPhaseChangeTwoPhaseMixture() {}
};

} // End namespace Foam

#endif
