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
    dimensionedScalar("eps",dimless,1.0)
);


Info<< "Reading diffusivity DT\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT("DT", dimViscosity, transportProperties);

bool hasTensorD = IOobject
(
    "D",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::NO_WRITE
).typeHeaderOk<volTensorField>();

autoPtr<volTensorField> DPtr;

if (hasTensorD)
{
    Info << "Reading tensorial diffusion field D\n" << endl;
    DPtr.reset(new volTensorField
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ));
}
else
{
    Info << "No tensor D found — using isotropic DT * I\n" << endl;
}

volVectorField J
(
     IOobject
     (
          "J",
          runTime.timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE
     ),
     -eps*(DT * fvc::grad(T))
);

if (hasTensorD)
{
	J = -eps*((*DPtr) & fvc::grad(T));
}

#include "createFvOptions.H"
