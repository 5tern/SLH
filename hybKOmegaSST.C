/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hybKOmegaSST.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSmooth.H"
#include "bound.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(hybKOmegaSST, 0);
addToRunTimeSelectionTable
(
    RAStransportModelIncompressibleTurbulenceModel,
    hybKOmegaSST,
    dictionary
);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

dimensionedScalar hybKOmegaSST::cD(const volSymmTensorField& D) const
{

    volSymmTensorField LL = dev(filter_(sqr(U())) - (sqr(filter_(U()))));

    volSymmTensorField MM =
      sqr(delta())*(filter_(mag(D)*(D)) - 4*mag(filter_(D))*filter_(D));

    dimensionedScalar MMMM = average(magSqr(MM));

    if (MMMM.value() > VSMALL)
    {
        return average(LL && MM)/MMMM;
    }
    else
    {
        return 0.0;
    }
}

tmp<volScalarField> hybKOmegaSST::F1(const volScalarField& CDkOmega) const
{
    volScalarField CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    volScalarField arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

tmp<volScalarField> hybKOmegaSST::F2() const
{
    volScalarField arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}

// both rd and fd should be const methods, but 
// for debug purposes I need them to be able
// to write fd, rd

tmp<volScalarField> hybKOmegaSST::rd
(
    const volScalarField& visc,
    const volScalarField& S,
    const volScalarField& Vorticity
)
{
    return
    visc
    /(
        max
        (
            sqrt(0.5*(S*S + Vorticity*Vorticity)),
            dimensionedScalar("SMALL", S.dimensions(), SMALL)
        )*sqr(0.41*y_)
      + dimensionedScalar
        (
            "ROOTVSMALL",
            dimensionSet(0, 2 , -1, 0, 0),
            1e-10
        )
    );

}


tmp<volScalarField> hybKOmegaSST::fd
(
    const volScalarField& visc,
    const volScalarField& S,
    const volScalarField& Vorticity
)
{
#ifdef debugModel
    rd_ = rd(visc, S, Vorticity);
    fd_ = 1 - tanh(pow(8.0*rd(visc, S, Vorticity), 3.0));
#endif
    return 1 - tanh(pow(8.0*rd(visc, S, Vorticity), 3.0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hybKOmegaSST::hybKOmegaSST
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type 
)
:
    RAStransportModelIncompressibleTurbulenceModel
    (
            type,
            alpha,
            rho,
            U,
            alphaRhoPhi,
            phi,
            transport,
            propertiesName
    ),
    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    
    Cint_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cint",
            coeffDict_,
            1.0
        )
    ),

    CN_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CN",
            coeffDict_,
            1.0
        )
    ),
       
    x1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "x1",
            coeffDict_,
            0.95
        )
    ),  
     
    x2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "x2",
            coeffDict_,
            1.05
        )
    ),

    dFilterIter_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "dFilterIter",
            coeffDict_,
            0
        )
    ),

 
    omegaSmall_
    (
        "omegaSmall", 
         dimensionSet(0, 0, -1, 0, 0), 
         SMALL
    ),
    
    shielding
    (
        Switch::lookupOrAddToDict
        (
            "shielding",
            coeffDict_,
            true
        )
    ),
    
    y_(wallDist::New(mesh_).y()),

    d_
    (
        IOobject
        (
            "d",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("d",dimLength,SMALL)
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
       mesh_
    ),

    rans_
    (
        IOobject
        (
            "rans",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("rans", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

#ifdef debugModel
    fd_
    (
        IOobject
        (
            "fd",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("fd", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rd_
    (
        IOobject
        (
            "rd",
            runTime_.timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("rd", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
#endif 

    simpleFilter_(this->mesh_),
    filterPtr_(LESfilter::New(this->mesh_, this->coeffDict())),
    filter_(filterPtr_()),
    delta_(LESdelta::New("delta", *this, this->coeffDict()))
    
{
    bound(omega_, omegaMin_);

    nut_ = a1_*k_/max(a1_*omega_, F2()*mag(symm(fvc::grad(U_))));
    nut_.correctBoundaryConditions();

    printCoeffs("hybKOmegaSST");
    
    const pointField& ps = mesh_.points();
    const edgeList& es = mesh_.edges();
    const labelListList& eCells = mesh_.cellEdges();

    scalar edgeLength = 0.0 ;
    volScalarField maxEdgeLength(d_);
    dimensionedScalar l("1m",dimLength,1.0);
    
    forAll( maxEdgeLength , I )
    {
        edgeLength = 0.0;
        forAll( eCells[I], J )
        {
	        if ( ( es[eCells[I][J]].mag( ps ) ) > edgeLength ) 
	            edgeLength = es[ eCells[I][J] ].mag( ps );
        }
        maxEdgeLength[I] = edgeLength; //*l;
    }
   
    forAll( d_ , I )
    {
        d_[I] = sqrt( 0.5*( maxEdgeLength[I]*maxEdgeLength[I] + delta()[I]*delta()[I] ) ); 
    }

    for(label i=0; i < static_cast<label>(dFilterIter_.value()) ; i++)
    {
        Info<<"d smoothing, iteration " << i << endl;
        d_ = filter_(d_);
    }
    d_.boundaryFieldRef() = y_.boundaryField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> hybKOmegaSST::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}

void hybKOmegaSST::correctNut
(
    const volScalarField& S2
)
{
    nut_ = a1_*k_/max(a1_*omega_, F2()*sqrt(S2));
    volScalarField lPrandtl = CN_ *((100.0/9.0)*(pow(k_,0.5)/(omega_ + omegaSmall_))); // integral length scale using Prandtl formula
    
    volSymmTensorField D = dev(symm(fvc::grad(U_)));
    volScalarField nuSGS = mag(cD(D))*sqr(rans_*d_+(1-rans_)*delta())*sqrt(2*magSqr(D)); // SGS stress using DSM

    if(shielding)
    {
        lPrandtl = lPrandtl * fd(nut_*rans_ + nuSGS*(scalar(1.0)-rans_) + nu(),sqrt(2.0)*
        mag(symm(fvc::grad(U_))),sqrt(2.0)*mag(skew(fvc::grad(U_))));
    }

    scalar h  = 0;
    
    scalar x1 = x1_.value(); 
    scalar x2 = x2_.value();

    // Coefficient for a tanh function, 
    // ensuring, that at x2, blending function will be equal 1-eps
    // the lower eps is, the steeper the blending function is
    
    scalar eps = 1e-12;
    
    // Coefficient for a tanh function
    scalar c= 2*(atanh(2*(1-eps)-1.0)/(x2-x1)); 
    
    forAll( d_, I )
    {
	    rans_[I] = 1; 
            h = lPrandtl[I] / d_[I] ;
	    h /= Cint_.value();    // how many cells to resolve a vortex?

	    if ( h  > x2 ) 
	    {
	        nut_[I] = nuSGS[I];   // cell is in LES regio
	        rans_[I] = 0;
	    }
	    else
	    {
	        if ( h > x1 ) // cell is buffer region between URANS and LES (see hybrid model desc.)
	        {
             
                     rans_[I]  =  0.5 - 0.5 * tanh(c*( h - (x1+x2)/2) );
                     nut_[I] = nut_[I]*rans_[I] + nuSGS[I]*(1.0-rans_[I]);
	        }
	    }        
    }

    nut_.correctBoundaryConditions();
}


void hybKOmegaSST::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}

tmp<volSymmTensorField> hybKOmegaSST::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> hybKOmegaSST::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}

tmp<fvVectorMatrix> hybKOmegaSST::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

bool hybKOmegaSST::read()
{
    if (RASModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        Cint_.readIfPresent(coeffDict());
        CN_.readIfPresent(coeffDict());
        x1_.readIfPresent(coeffDict());
        x2_.readIfPresent(coeffDict());
        
        return true;
    }
    else
    {
        return false;
    }
}


void hybKOmegaSST::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }
    
    volTensorField gradU = fvc::grad(U_);
    
    const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    volScalarField G(GName(), nut_*S2);

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega =
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_;

    volScalarField F1 = this->F1(CDkOmega);

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)*S2
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );
    
    omegaEqn.ref().relax();
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    bound(omega_, this->omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G, c1_*betaStar_*k_*omega_)
      - fvm::Sp(betaStar_*omega_, k_)
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, this->kMin_);

	correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
