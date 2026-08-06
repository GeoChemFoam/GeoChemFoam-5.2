{    
    R = 0*R;
    forAll(kineticPhaseReactions,j)
    {
        dimensionedScalar Mws = speciesMixture.Mws(j);
        dimensionedScalar rhos = speciesMixture.rhos(j);
        volScalarField Rj = speciesMixture.Rate(j);
        volScalarField wsj = 1-eps;
        forAll(kineticPhases,k)
        {
            if (kineticPhaseReactions[j]==kineticPhases[k])
            {
                volScalarField ws    = speciesMixture.Ws(k);
                volScalarField omega = speciesMixture.omega(k);
                Rj *= (1-omega);
                wsj*=ws;
            }
        }

        surfaceScalarField wsjf = fvc::interpolate(wsj)/(fvc::interpolate(1-eps+1e-15));
        surfaceScalarField phiR = Mws/rhos*fvc::interpolate(Rj)*wsjf*nEpsf;

        surfaceScalarField phiREps = fvc::flux(phiR,eps,"div(phiR,eps)");

        R+=eps*fvc::div(phiR)-fvc::div(phiREps);
    }


    forAll(solutionSpecies, i)
    {
	volScalarField& Yi = speciesMixture.Y(i);
        dimensionedScalar DYi = speciesMixture.DY(i);
        bool hasTensorDYi = speciesMixture.hasTensorDY(i);
        Ri = 0*Ri;

        forAll(kineticPhaseReactions,j)
        {
            volScalarField Rj = speciesMixture.Rate(j);
            const scalarList& scoeffj = speciesMixture.scoeff(j);
        
            forAll(kineticPhases,k)
            {
                if (kineticPhaseReactions[j]==kineticPhases[k])
                {
                    volScalarField omega = speciesMixture.omega(k);
                    Rj *= (1-omega);
                }
            }

            surfaceScalarField phiR = -scoeffj[i]*fvc::interpolate(Rj)*nEpsf;

            surfaceScalarField phiREps = scoeffj[i]*fvc::interpolate(Rj)*fvc::flux(-nEpsf,eps,"div(nEpsf,eps)");
            Ri += fvc::div(phiREps)-eps*fvc::div(phiR);
        }

        volScalarField D = eps*DYi;
                
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix YiEqn
            (
	        fvm::ddt(eps,Yi)
                + fvm::div(phi, Yi, "div(phi,Yi)")
                + fvm::SuSp(-fvc::div(phi),Yi)
                + Ri
            );

            if (hasTensorDYi)
            {
                const volTensorField& DTensorYi = speciesMixture.DTensorY(i);
                YiEqn += - fvm::laplacian(eps*DTensorYi, Yi);
            }
            else
            {
                YiEqn += - fvm::laplacian(eps*DYi, Yi);
            }

            YiEqn.relax();
            YiEqn.solve(mesh.solver("Yi"));
            Yi = max(Yi,0*Yi);
        }

	Info<< solutionSpecies[i] << " concentration = "
			<< "  Min(Yi) = " << gMin(Yi)
			<< "  Max(Yi) = " << gMax(Yi)
			<< endl;
    }
}

speciesMixture.correct();
