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
    processConcentration 

Description
    calculate average concentration of each species for each time-step

\*---------------------------------------------------------------------------*/


#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "fvCFD.H"
#include "argList.H"
#include "primitivePatchInterpolation.H"
#include "timeSelector.H"
#include "speciesTable.H"
#include "upwind.H"


using namespace Foam;


//#include "timeDataList.H"


int main(int argc, char *argv[])
{
	#   include "setRootCase.H"
	#   include "createTime.H"
	#   include "createNamedMesh.H"


	instantList timeList = timeSelector::select0(runTime, args);
//	timeDataList dataList(timeList.size());


	IOdictionary postProcessDict 
	(
		IOobject
		(
			"postProcessDict",
			"system",
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::NO_WRITE
		)
	);

	scalar z1
	(
		postProcessDict.lookupOrDefault<scalar>("z1",-1e30)
	);


        scalar z2
        (
                postProcessDict.lookupOrDefault<scalar>("z2",1e30)
        );

        scalar y1
        (
                postProcessDict.lookupOrDefault<scalar>("y1",-1e30)
        );


        scalar y2
        (
                postProcessDict.lookupOrDefault<scalar>("y2",1e30)
        );


        scalar x1
        (
                postProcessDict.lookupOrDefault<scalar>("x1",-1e30)
        );


        scalar x2
        (
                postProcessDict.lookupOrDefault<scalar>("x2",1e30)
        );

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
        IOdictionary thermoPhysicalProperties
        (
                IOobject
                (
                       "thermoPhysicalProperties",
                       "constant",
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                )
        );

	const dictionary& solutionSpeciesDict = thermoPhysicalProperties.subDict("solutionSpecies");
	speciesTable solutionSpecies
	(
		solutionSpeciesDict.toc()
	);


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
             dimensionedScalar("eps0",dimless,1.0)
        );

        volScalarField aI
        (
            IOobject
            (
                "aI",
                 runTime.timeName(),
                 mesh,
                 IOobject::READ_IF_PRESENT,
                 IOobject::NO_WRITE
             ),
             mesh,
             dimensionedScalar("aI",dimless/dimLength,0.0)
        );

        dimensionedScalar deltaN = 1e-13/pow(average(mesh.V()), 1.0/3.0);

        volVectorField gradEps = fvc::grad(eps0);
        surfaceVectorField gradEpsf = fvc::interpolate(gradEps);
        surfaceVectorField nEpsv = -gradEpsf/(mag(gradEpsf) + deltaN);
        surfaceScalarField  nEpsf = nEpsv & mesh.Sf();

        surfaceScalarField eps_f = upwind<scalar>(mesh, -nEpsf).interpolate(eps0);

        volScalarField aSurf0
        (
            IOobject
            (
                "aSurf",
                 runTime.timeName(),
                 mesh,
                 IOobject::READ_IF_PRESENT,
                 IOobject::NO_WRITE
             ),
             aI-fvc::div(eps_f*nEpsf)+eps0*fvc::div(nEpsf)
        );


        forAll(solutionSpecies,i)
        {
            string filename = solutionSpecies[i]+"_Conc.csv";
	    std::ofstream csvfile(filename.c_str());
	    csvfile << "time C\n";

        
	    forAll(timeList, timeStep)
	    {	
	    	runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;

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


              	PtrList<volScalarField> Y(solutionSpecies.size());

		Y.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					solutionSpecies[i],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);


		volScalarField& Yi = Y[i];

                scalarField& vol = clip;

		scalar T = gSum(eps*Yi*vol)/(gSum(eps*vol)+1e-30);

		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << T << "\n";
		}
	    }

	    csvfile.close();
        }

        if (thermoPhysicalProperties.found("surfaceSpecies"))
        {
	    const dictionary& surfaceSpeciesDict = thermoPhysicalProperties.subDict("surfaceSpecies");
	    speciesTable surfaceSpecies
	    (
		surfaceSpeciesDict.toc()
	    );

            forAll(surfaceSpecies,i)
            {
            string filename = surfaceSpecies[i]+"_Conc.csv";
	    std::ofstream csvfile(filename.c_str());
	    csvfile << "time S\n";

        
	    forAll(timeList, timeStep)
	    {	
	    	runTime.setTime(timeList[timeStep], timeStep);

		Info<< endl<<timeStep<< "    Time = " << runTime.timeName() << endl;

                volScalarField aSurf
                (
                        IOobject
                        (
                                "eps",
                                runTime.timeName(),
                                mesh,
                                IOobject::READ_IF_PRESENT,
                                IOobject::NO_WRITE
                        ),
                        aSurf0
                );


              	PtrList<volScalarField> sY(surfaceSpecies.size());

		sY.set
		(
			i,
			new volScalarField
			(
				IOobject
				(
					surfaceSpecies[i],
					mesh.time().timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE
				),
				mesh
			)
		);


		volScalarField& Yi = sY[i];

                scalarField& vol = clip;

		scalar sumT = gSum(aSurf*Yi*vol);

                scalar sumA = gSum(aSurf*vol);

                forAll(Yi.boundaryFieldRef(), patchi)
                {
                    if (Yi.boundaryFieldRef()[patchi].type()=="reactingWall")
                    {
                        const surfaceScalarField& magSf = mesh.magSf();
                        scalarField& Yfaces   = Yi.boundaryFieldRef()[patchi];
                        sumT =sumT+gSum(Yfaces*magSf.boundaryField()[patchi]);
                        sumA =sumA+gSum(magSf.boundaryField()[patchi]);
                    }
                }

                scalar sT = sumT/(sumA+1e-30);
		if (Pstream::master())
		{
			
			csvfile << runTime.timeName() << " " << sT << "\n";
		}
	    }

	    csvfile.close();
        }

        }

	return 0;
}


// ************************************************************************* //
