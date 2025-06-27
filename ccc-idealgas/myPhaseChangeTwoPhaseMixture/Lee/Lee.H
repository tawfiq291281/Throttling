#ifndef Lee_H
#define Lee_H

#include "myPhaseChangeTwoPhaseMixture.H"
#include "volFields.H"
#include "dimensionedScalar.H"

namespace Foam
{
    class Lee
    :
        public myPhaseChangeTwoPhaseMixture
    {
    protected:
        dimensionedScalar r_;
        dimensionedScalar L_;
        dimensionedScalar a_;
        dimensionedScalar b_;
        dimensionedScalar c_;
        volScalarField T_;

    public:
        TypeName("Lee");

        Lee
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            const volScalarField& alpha1,
            const volScalarField& alpha2
        );

        virtual ~Lee() {}

        virtual const volScalarField& T() const;

        virtual tmp<volScalarField> heSource() const;

        virtual tmp<volScalarField> T(const volScalarField& h) const;

        virtual tmp<volScalarField> mDotAlphal() const;

        virtual tmp<volScalarField> vDotP() const;

        virtual tmp<volScalarField> rho2() const;

        virtual void correct();

        virtual bool read();

        // Accessors for Antoine coefficients
        inline const dimensionedScalar& getA() const { return a_; }
        inline const dimensionedScalar& getB() const { return b_; }
        inline const dimensionedScalar& getC() const { return c_; }

        // New accessors for phase change rates
        virtual tmp<volScalarField> mDotEvaporation() const;
        virtual tmp<volScalarField> mDotCondensation() const;
    };
}

#endif
