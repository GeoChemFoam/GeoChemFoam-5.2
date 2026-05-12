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
    processFF

Description
    calculate porosity and formation factor

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "fvCFD.H"
#include "argList.H"
#include "primitivePatchInterpolation.H"
#include "timeSelector.H"
#include "tractionDisplacementFvPatchVectorField.H"


using namespace Foam;

int main(int argc, char *argv[])
{
	#   include "setRootCase.H"
	#   include "createTime.H"
        #   include "createNamedMesh.H"

	instantList timeList = timeSelector::select0(runTime, args);
 
        scalar poro = 0.0;
	scalar Lx = 0.0;
	scalar Ly = 0.0;
	scalar Lz = 0.0;
        scalar stressx = 0.0;
        scalar stressy =0.0;
        scalar stressz = 0.0;
        scalar strainx =0.0;
        scalar strainy = 0.0;
        scalar strainz = 0.0; 

	std::ofstream csvfile("StressStrain.csv");
	csvfile << "time poro stressx stressy stressz strainx strainy strainz\n";

        volScalarField eps0
        (
            IOobject
            (
                "eps",
                 runTime.timeName(),
                 mesh,
                 IOobject::READ_IF_PRESENT,
                 IOobject::NO_WRITE
             ),
             mesh,
             dimensionedScalar("eps",dimless,0.0)
        );

        IOdictionary postProcessDict
        (
                IOobject
                (
                        "postProcessDict",
                        "system",
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                )
        );

        scalar z1
        (
            readScalar(postProcessDict.lookup("z1"))
        );


        scalar z2
        (
            readScalar(postProcessDict.lookup("z2"))
        );

        scalar y1
        (
            readScalar(postProcessDict.lookup("y1"))
        );


        scalar y2
        (
            readScalar(postProcessDict.lookup("y2"))
        );


        scalar x1
        (
            readScalar(postProcessDict.lookup("x1"))
        );


        scalar x2
        (
            readScalar(postProcessDict.lookup("x2"))
        );

        word leftName
        (
            postProcessDict.lookup("leftName")
        );

        label leftPatchID  = mesh.boundaryMesh().findPatchID(leftName);

        word rightName
        (
            postProcessDict.lookup("rightName")
        );

        label rightPatchID  = mesh.boundaryMesh().findPatchID(rightName);

        word bottomName
        (
            postProcessDict.lookup("bottomName")
        );

        label bottomPatchID  = mesh.boundaryMesh().findPatchID(bottomName);

        word topName
        (
            postProcessDict.lookup("topName")
        );

        label topPatchID  = mesh.boundaryMesh().findPatchID(topName);

        word frontName
        (
            postProcessDict.lookup("frontName")
        );

        label frontPatchID  = mesh.boundaryMesh().findPatchID(frontName);


        word backName
        (
            postProcessDict.lookup("backName")
        );

        label backPatchID  = mesh.boundaryMesh().findPatchID(backName);


        Info<< "Reading mechanical properties\n" << endl;

        IOdictionary mechanicalProperties
        (
            IOobject
            (
                "mechanicalProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

        const dictionary& rhoDict(mechanicalProperties.subDict("rho"));
        word rhoType(rhoDict.get<word>("type"));

        autoPtr<volScalarField> rhoPtr;

        IOobject rhoIO
        (
            "rho",
            runTime.timeName(0),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        );

        if (rhoType == "uniform")
        {
            scalar rhoValue(rhoDict.get<scalar>("value"));

            rhoPtr.reset
            (
                new volScalarField
                (
                    rhoIO,
                    mesh,
                    dimensionedScalar
                    (
                        "rho",
                        dimMass/dimVolume,
                        rhoValue
                    )
                )
            );
        }
        else if (rhoType == "field")
        {
            rhoIO.readOpt(IOobject::MUST_READ);

            rhoPtr.reset
            (
                new volScalarField
                (
                    rhoIO,
                    mesh
                )
            );
        }
        else
        {
           FatalErrorInFunction
                << "Valid type entries are uniform or field for rho"
                << abort(FatalError);
        }

        volScalarField& rho = rhoPtr();


	forAll(timeList, timeStep)
	{
		runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;
	
		#   include "createNamedMesh.H"

                volScalarField clip
                (
                   IOobject
                   (
                        "clip",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimensionedScalar("clip",dimVolume,0.0)
                );


                volScalarField coordz=mesh.C().component(2);
                volScalarField coordy=mesh.C().component(1);
                volScalarField coordx=mesh.C().component(0);

                forAll(mesh.cells(),j)
                {
                    scalar zj=coordz[j];
                    scalar yj=coordy[j];
                    scalar xj=coordx[j];
                    if (zj>=z1 && zj<=z2 && yj>=y1 && yj<=y2 && xj>=x1 && xj<=x2)
                    {
                            clip[j]=mesh.V()[j];
                    }
                }

                volScalarField eps
                (
                        IOobject
                        (
                                "eps",
                                runTime.timeName(),
                                mesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::NO_WRITE
                        ),
                        eps0
                );

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

                volScalarField epsSolid = 1-eps;
                poro = (1.0-epsSolid.weightedAverage(clip).value()*gSum(clip)/(x2-x1)/(y2-y1)/(z2-z1));

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
                    mesh,
                    dimensionedSymmTensor("sigmaD",dimPressure/dimDensity,symmTensor::zero)
                );

                volSymmTensorField sigma = rho*sigmaD;

                volScalarField sigmaxx = sigma.component(symmTensor::XX);
                volScalarField sigmayy = sigma.component(symmTensor::YY);
                volScalarField sigmazz = sigma.component(symmTensor::ZZ);

                volScalarField Dx = D.component(0);
                volScalarField Dy = D.component(1);
                volScalarField Dz = D.component(2);

                scalarField vol = mesh.V();

                const fvPatchScalarField& DRight = Dx.boundaryField()[rightPatchID];
                const fvPatchScalarField& DLeft  = Dx.boundaryField()[leftPatchID];
                scalarField magSfRight = DRight.patch().magSf();
                scalarField magSfLeft  = DLeft.patch().magSf();
		const vectorField& CfR = mesh.Cf().boundaryField()[rightPatchID];
                const vectorField& CfL = mesh.Cf().boundaryField()[leftPatchID];
		Lx=gMax(CfR.component(0))-gMin(CfL.component(0));

                const fvPatchScalarField& DTop     = Dy.boundaryField()[topPatchID];
                const fvPatchScalarField& DBottom  = Dy.boundaryField()[bottomPatchID];
                scalarField magSfTop     = DTop.patch().magSf();
                scalarField magSfBottom  = DBottom.patch().magSf();
                const vectorField& CfT = mesh.Cf().boundaryField()[topPatchID];
                const vectorField& CfB = mesh.Cf().boundaryField()[bottomPatchID];
                Ly=gMax(CfT.component(1))-gMin(CfB.component(1));

                const fvPatchScalarField& DBack = Dz.boundaryField()[backPatchID];
                const fvPatchScalarField& DFront  = Dz.boundaryField()[frontPatchID];
                scalarField magSfBack   = DBack.patch().magSf();
                scalarField magSfFront  = DFront.patch().magSf();
		const vectorField& CfD = mesh.Cf().boundaryField()[backPatchID];
                const vectorField& CfF = mesh.Cf().boundaryField()[frontPatchID];
                Lz=gMax(CfD.component(2))-gMin(CfF.component(2));

                stressx=gSum(sigmaxx*vol)/Lx/Ly/Lz;
                stressy=gSum(sigmayy*vol)/Lx/Ly/Lz;
                stressz=gSum(sigmazz*vol)/Lx/Ly/Lz;


		strainx = (gSum(DRight*magSfRight)/gSum(magSfRight)-gSum(DLeft*magSfLeft)/gSum(magSfLeft))/Lx;
		strainy = (gSum(DTop*magSfTop)/gSum(magSfTop)-gSum(DBottom*magSfBottom)/gSum(magSfBottom))/Ly;
                strainz = (gSum(DBack*magSfBack)/gSum(magSfBack)-gSum(DFront*magSfFront)/gSum(magSfFront))/Lz;

		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << poro << " " << stressx << " " << stressy << " " << stressz << " " << strainx << " " << strainy << " " << strainz << " " <<"\n";
		}
	}

	csvfile.close();

	return 0;
}


// ************************************************************************* //
