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

#include "createPhi.H"

volVectorField B
(
    IOobject
    (
        "B",
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


Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT("DT", transportProperties);

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

dimensionedVector Umoy ("Umoy", U.weightedAverage(mesh.V()));

dimensionedScalar poro ("poro",eps.weightedAverage(mesh.V()));

Info<< "Umoy = " << Umoy.value() << nl << endl;

Info << "poro = " << poro.value() << nl << endl;


#include "createFvOptions.H"


