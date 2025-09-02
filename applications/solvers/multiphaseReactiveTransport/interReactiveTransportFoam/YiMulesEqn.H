{
    word alpharScheme("div(phirb,alpha)");
    word YiScheme("div(phi,Yi)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(pSpeciesMixture().cYi()*mag(phi/mesh.magSf()));

    surfaceScalarField phir(phic*mixture.nHatf());

    basicTwoPhaseMultiComponentTransportMixture& speciesMixture =pSpeciesMixture();
    const speciesTable& solutionSpecies = speciesMixture.species();

    surfaceScalarField phiComp = fvc::flux
    (
        -fvc::flux(-phir, alpha2, alpharScheme),
        alpha1,
        alpharScheme
    );

    forAll(solutionSpecies, i)
    {
        volScalarField& Yi = speciesMixture.Y(i);

        //compute diffusion coefficient
        surfaceScalarField D = speciesMixture.DmY(i);
        
        //calculate Henry's transfer flux
        surfaceScalarField phiH = speciesMixture.phiH(i);
        surfaceScalarField phiHUp = speciesMixture.phiHUp(i);
        surfaceScalarField phiHDown = speciesMixture.phiHDown(i);

        scalar maxYi = max(gMax(Yi.oldTime()),gMax(Yi.oldTime().boundaryField()))+1e-30;

        Yi.oldTime() == Yi.oldTime()/maxYi;
        Yi == Yi/maxYi;
         
        surfaceScalarField phiDY = -D*fvc::snGrad(Yi)*mesh.magSf();

        if (phiHScheme == "Gauss linear") phiDY += fvc::flux(phiH,Yi,"div(phiH,Yi)");
        else if (phiHScheme == "Gauss upwind")
        {
            phiDY += fvc::flux(phiHUp,Yi,"div(phiH,Yi)")
            + fvc::flux(phiHDown,Yi,"div(phiH,Yi)");
        }
        else
        {
            Info<< "div(phiH,Yi) should be equal to Gauss linear or Gauss upwind"
                            << endl
                            << abort(FatalError);
        }

        volScalarField SpY = -fvc::div(phiDY);

        tmp<surfaceScalarField> tYiPhi1Un
        (
            fvc::flux
            (
                phi,
                Yi,
                YiScheme
            )
            +   phiComp*speciesMixture.compressionCoeff(i)
        );

        {
            surfaceScalarField YiPhi10 = tYiPhi1Un;

            MULES::explicitSolve
            (
                geometricOneField(),
                Yi,
                phi,
                YiPhi10,
                zeroField(),
                SpY,
                oneField(),
                zeroField()
            );
        }

        Yi.oldTime() == Yi.oldTime()*maxYi;
        Yi == Yi*maxYi;

        Info<< "Species concentration = "
        << Yi.weightedAverage(mesh.V()).value()
        << "  Min(Yi) = " << gMin(Yi.internalField())
        << "  Max(Yi) = " << gMax(Yi.internalField())
        << endl;

    }
 }
