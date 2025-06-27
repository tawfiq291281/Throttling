#include "myPhaseChangeTwoPhaseMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"

namespace Foam
{
    defineTypeNameAndDebug(myPhaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(myPhaseChangeTwoPhaseMixture, components);

    myPhaseChangeTwoPhaseMixture::myPhaseChangeTwoPhaseMixture
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& alpha1,
        const volScalarField& alpha2
    )
    :
        incompressibleTwoPhaseMixture(U, phi),
        L_
        (
            U.db().lookupObject<IOdictionary>("transportProperties").get<dimensionedScalar>("L")
        ),
        cpl1_
        (
            U.db().lookupObject<IOdictionary>("transportProperties").subDict("liquid").get<dimensionedScalar>("cp")
        ),
        cpv2_
        (
            U.db().lookupObject<IOdictionary>("transportProperties").subDict("vapour").get<dimensionedScalar>("cp")
        ),
        alpha1_(alpha1),
        alpha2_(alpha2),
        phaseChangeDict_
        (
            U.db().lookupObject<IOdictionary>("transportProperties").subDict("phaseChangeDict")
        )
    {}

    autoPtr<myPhaseChangeTwoPhaseMixture> myPhaseChangeTwoPhaseMixture::New
    (
        const volVectorField& U,
        const surfaceScalarField& phi,
        const volScalarField& alpha1,
        const volScalarField& alpha2
    )
    {
        word phaseChangeTwoPhaseMixtureType
        (
            U.db().lookupObject<IOdictionary>("transportProperties").get<word>("phaseChangeTwoPhaseMixture")
        );

        auto cstrIter = componentsConstructorTablePtr_->find(phaseChangeTwoPhaseMixtureType);

        if (!cstrIter.found())
        {
            FatalErrorInFunction
                << "Unknown phaseChangeTwoPhaseMixtureType type "
                << phaseChangeTwoPhaseMixtureType << endl
                << "Valid phaseChangeTwoPhaseMixture types are : " << endl
                << componentsConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        return cstrIter()(U, phi, alpha1, alpha2);
    }
}
