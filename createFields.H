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


    Info<< "Reading transportProperties\n" << endl;

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

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );


    Info<< "Reading materialProperties\n" << endl;

    IOdictionary materialProperties
    (
        IOobject
        (
            "materialProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading thermal conductivity\n" << endl;

    dimensionedScalar k
    (
        materialProperties.lookup("k")
    );
    
    
    volScalarField sourceTerm
    (
        IOobject
        (
            "sourceTerm",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("sourceTerm", dimEnergy/dimTime/dimVolume, 0.0)  //source term value is set in the solver myLaplacianFoam.C
    );

