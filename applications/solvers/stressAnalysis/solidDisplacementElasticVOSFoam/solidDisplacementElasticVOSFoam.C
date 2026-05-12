/*---------------------------------------------------------------------------*\

License
    This file is part of GeoChemFoam, an Open source software using OpenFOAM
    for multiphase multicomponent reactive transport simulation in pore-scale
    geological domain.

    GeoChemFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version. See <http://www.gnu.org/licenses/>.

    The code was developed by Dr Julien Maes as part of his research work for
    the GeoChemFoam Group at Heriot-Watt University. Please visit our
    website for more information <https://github.com/GeoChemFoam>.

Application
    solidDisplacementElasticVOSFoam

Group
    grpStressAnalysisSolvers

Description
    Transient segregated finite-volume solver of linear-elastic,
    small-strain deformation of a solid body, with optional thermal
    diffusion and thermal stresses.

    Simple linear elasticity structural analysis code.
    Solves for the displacement vector field D, also generating the
    stress tensor field sigma.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "upwind.H"
#include "downwind.H"
#include "tractionDisplacementFvPatchVectorField.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient segregated finite-volume solver of linear-elastic,"
        " small-strain deformation of a solid body"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControls.H"

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating displacement field\n" << endl;

    while (simple.loop())
    {
        Info<< "Time: " << runTime.value() << nl << endl;

        #include "readSolidDisplacementFoamControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;

        do
        {
            {
                fvVectorMatrix DEqn
                (
                    fvm::d2dt2(D)
                 ==
                    fvm::laplacian(2*mu+lambda,D, "laplacian(DD,D)")
                  //immersed boundary correction
		  - fvm::div(fvc::interpolate(2*mu+lambda)*phiEpsS,D,"div(sigmaD)")//zeroGradient correction on immersed boundary
		  + epsSolid*fvc::div(fvc::interpolate(2*mu+lambda)*mesh.magSf()*snGradDcorr)//snGradDcorr correction on immersed boundary 
                  +divSigmaExp
                );

                initialResidual = DEqn.solve().max().initialResidual();
                D.relax();
            }

            {
                
                gradD=fvc::grad(D);

                //Displacement interpolated on face + immersed boundary correction
                DS = fvc::interpolate(D)/fvc::interpolate(epsSolid+1e-30)+snGradDcorr*deltaXS;

                //gradient displacement at cell center
                gradDS = epsSolid*fvc::grad(DS);

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
		snGradDS = fvc::snGrad(D) - fvc::snGrad(epsSolid)*fvc::interpolate(D)/fvc::interpolate(epsSolid+1e-30)+snGradDcorr;

                //gradient displacement at face center
                gradDSf = fvc::interpolate(epsSolid*gradDS)/fvc::interpolate(epsSolid+1e-30);

                //correction to keep grad=snGrad in the direction of Sf
                gradDSf = gradDSf + nf*(snGradDS - (nf&gradDSf));
                
                //stress at cell center
                sigmaD = mu*twoSymm(gradDS) + (lambda*I)*tr(gradDS);

                //stress at face center
		sigmaDf = fvc::interpolate(mu)*twoSymm(gradDSf) + fvc::interpolate(lambda)*(I*tr(gradDSf));

                //forces on face
                Fsigma = mesh.Sf() & sigmaDf;

		#include "correctStressBoundary.H"

                //explicit divergence term
                divSigmaExp = fvc::surfaceIntegrate(Fsigma)
                            - fvc::laplacian((2*mu+lambda),D,"laplacian(DD,D)")
                            //immersed boundary correction
                            - Psigma //fluid pressure
                            + fvc::div(fvc::interpolate(2*mu+lambda)*phiEpsS,D,"div(sigmaD)") //zero gradient correction on immersed boundary
                            - epsSolid*fvc::div(fvc::interpolate(2*mu+lambda)*mesh.magSf()*snGradDcorr); //snGradCorr on immersed boundary

                            
                snGradDcorr = snGradDcorr -deltaI*((nf & sigmaDf) + fvc::interpolate(p)/fvc::interpolate(rho)/fvc::interpolate(eps+1e-30)*nf)/fvc::interpolate(2*mu+lambda);

            }

        } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);

        //#include "calculateStress.H"

        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
