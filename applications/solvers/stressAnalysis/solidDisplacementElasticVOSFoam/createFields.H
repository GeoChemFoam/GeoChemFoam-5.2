#include "readMechanicalProperties.H"

Info<< "Reading field D\n" << endl;
volVectorField D
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
);

//local porosity
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
    dimensionedScalar("eps",dimless,0.0)
);

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    pf*eps
);


Info<< "Reading field snGradDcorr\n" << endl;
surfaceVectorField snGradDcorr 
(
    IOobject
    (
      "snGradDcorr",
      runTime.timeName(),
      mesh,
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("snGradDcorr",dimless,Zero)
);

snGradDcorr.setOriented(true);

//volume fraction of solid
volScalarField epsSolid = 1 - eps;

//Flux of gradient epsSolid - for immersed boundary term
surfaceScalarField phiEpsS(1/fvc::interpolate(epsSolid+1e-30)*fvc::snGrad(epsSolid)*mesh.magSf());

//fluid pressure stress
volVectorField Psigma = fvc::reconstruct(fvc::interpolate(p)/fvc::interpolate(rho)/fvc::interpolate(eps+1e-30)*fvc::snGrad(eps)*mesh.magSf());

//Interface localisation (upwind to keep it in fluid cells)
surfaceScalarField fluxDirSolid = fvc::snGrad(epsSolid)*mesh.magSf();
surfaceScalarField epsSolidUp = upwind<scalar>(mesh,fluxDirSolid).interpolate(epsSolid);
surfaceScalarField deltaI = (1-epsSolidUp/(epsSolidUp+1e-30))*fvc::interpolate(epsSolid)/fvc::interpolate(epsSolid+1e-30);

volTensorField gradD(fvc::grad(D));

//Displacement interpolated on face + immersed boundary correction
surfaceScalarField deltaXS = sign(fvc::snGrad(eps))*mag(mesh.delta())/2;
surfaceVectorField DS = fvc::interpolate(D)/fvc::interpolate(epsSolid+1e-30)+snGradDcorr*deltaXS;
//gradient displacement at cell center
volTensorField gradDS(epsSolid*fvc::grad(DS));


//In the case of traction displacement, the gradient is given on the boundary, and we should use that one
forAll(mesh.boundary(), patchi)
{
   // Check if the patch field type of D matches tractionDisplacement
   if (isA<tractionDisplacementFvPatchVectorField>(D.boundaryField()[patchi]))
   {
        gradDS.boundaryFieldRef()[patchi] = gradD.boundaryFieldRef()[patchi];
   }
}

//surface gradient of displacement with immersed boundary correction
surfaceVectorField snGradDS = fvc::snGrad(D) - fvc::snGrad(epsSolid)*fvc::interpolate(D)/fvc::interpolate(epsSolid+1e-30)+snGradDcorr;

//gradient displacement at face center
surfaceTensorField gradDSf = fvc::interpolate(epsSolid*gradDS)/fvc::interpolate(epsSolid+1e-30);

//correction to keep grad=snGrad in the direction of Sf
surfaceVectorField nf = mesh.Sf()/mesh.magSf();
gradDSf = gradDSf + nf*(snGradDS - (nf&gradDSf));

Info<< "Calculating stress field sigmaD\n" << endl;
volSymmTensorField sigmaD
(
    IOobject
    (
        "sigmaD",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mu*twoSymm(gradDS) + lambda*(I*tr(gradDS))
);


surfaceSymmTensorField sigmaDf
(
    IOobject
    (
        "sigmaDf",        
        runTime.timeName(),
        mesh,              
        IOobject::READ_IF_PRESENT, 
        IOobject::AUTO_WRITE 
    ),
    fvc::interpolate(mu)*twoSymm(gradDSf) + fvc::interpolate(lambda)*(I*tr(gradDSf))
);


//Stress flux
surfaceVectorField Fsigma = mesh.Sf() & sigmaDf;

#include "correctStressBoundary.H"

Info<< "Calculating explicit part of div(sigma) divSigmaExp\n" << endl;
volVectorField divSigmaExp
(
    IOobject
    (
        "divSigmaExp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::surfaceIntegrate(Fsigma)
    -fvc::laplacian((2*mu + lambda),D,"laplacian(D,DD)")
    //immersed boundary correction
    - Psigma //fluid pressure
    +fvc::div(fvc::interpolate(2*mu+lambda)*phiEpsS,D,"div(phiEps,D)")//zero gradient correction on immersed boundary
    - epsSolid*fvc::div(fvc::interpolate(2*mu+lambda)*mesh.magSf()*snGradDcorr)//snGradcorr on immersed boundary
);

//update snGradCorr according to immersed boundary correction
snGradDcorr = snGradDcorr -deltaI*((nf & sigmaDf) + fvc::interpolate(p)/fvc::interpolate(rho)/fvc::interpolate(eps+1e-30)*nf)/fvc::interpolate(2*mu+lambda);

mesh.setFluxRequired(D.name());
