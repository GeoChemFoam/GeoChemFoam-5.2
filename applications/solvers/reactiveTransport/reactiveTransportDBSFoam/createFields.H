Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
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

volScalarField eps
(
    IOobject
    (
        "eps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("eps",dimless,1)
);

volScalarField epsSmooth
(
    IOobject
    (
        "epsSmooth",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    eps
);

scalar cS = simple.dict().lookupOrDefault<scalar>("cS", 0.0);
label  nS = simple.dict().lookupOrDefault<label>("nS", 0);

for(int i=0;i<nS;i++)
{
     epsSmooth = cS*fvc::average(fvc::interpolate(epsSmooth))+(1-cS)*epsSmooth;
}

dimensionedScalar deltaN = 1e-13/pow(average(mesh.V()), 1.0/3.0);

volVectorField gradEps = fvc::grad(eps);
surfaceVectorField gradEpsf = fvc::interpolate(gradEps);
surfaceVectorField nEpsv = -gradEpsf/(mag(gradEpsf) + deltaN);
surfaceScalarField  nEpsf = nEpsv & mesh.Sf();

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

#include "createPhi.H"

Info << "Create species mixture\n" << endl;
autoPtr<basicReactiveTransportMixture> pSpeciesMixture
(
    basicReactiveTransportMixture::New(mesh)
);

basicReactiveTransportMixture& speciesMixture = pSpeciesMixture();
speciesMixture.correct();

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);


dimensionedScalar nu("nu",dimViscosity,transportProperties);


autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

dimensionedScalar kf(dimensionedScalar::getOrDefault("kf",transportProperties,dimless/sqr(dimLength),Zero));


volScalarField Kinv
(
    IOobject
    (
        "Kinv",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    kf*pow(1-eps,2)/pow(eps,3)
);


Info << "Reading reaction rate if present\n" << endl;
volScalarField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("R",dimless/dimTime, 0.0)
);

volScalarField Ri
(
    IOobject
    (
        "Ri",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMoles/dimVolume/dimTime,
        0.0
    )
);

#include "createMRF.H"
#include "createFvOptions.H"
