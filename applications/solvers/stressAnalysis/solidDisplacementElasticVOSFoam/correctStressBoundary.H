const volVectorField::Boundary& D_bound  = D.boundaryField();
const surfaceVectorField::Boundary& nf_bound  = nf.boundaryField();
const surfaceScalarField::Boundary& magSf_bound  = mesh.magSf().boundaryField();
const volScalarField::Boundary& rho_bound  = rho.boundaryField();
surfaceVectorField::Boundary& F_bound = Fsigma.boundaryFieldRef();

const fvBoundaryMesh& boundary = mesh.boundary();

forAll(boundary, patchi)
{
    if (isA<tractionDisplacementFvPatchVectorField>(D_bound[patchi]))
    {
        tractionDisplacementFvPatchVectorField& D_bound_cast =
        const_cast<tractionDisplacementFvPatchVectorField&>
        (
            refCast<const tractionDisplacementFvPatchVectorField>
            (
                D_bound[patchi]
            )
         );

         fvsPatchVectorField& F_bound_p = F_bound[patchi];
         const vectorField& traction = D_bound_cast.traction();
         const scalarField& pressure = D_bound_cast.pressure();
         const vectorField& nfbound = nf_bound[patchi];
         const scalarField& magSfbound = magSf_bound[patchi];
         const scalarField& rhobound = rho_bound[patchi];
          

         F_bound_p = magSfbound*(traction - pressure*nfbound)/rhobound;
    }
}

//correct force on immersed boundary
Fsigma = (1-deltaI)*Fsigma - deltaI*fvc::interpolate(epsSolid)*fvc::interpolate(p)/fvc::interpolate(rho)/fvc::interpolate(eps+1e-30)*nf*mesh.magSf();

