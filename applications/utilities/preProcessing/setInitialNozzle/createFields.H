Info<< "Reading thermophysical properties\n" << endl;

IOdictionary initialConditions
(
    IOobject
    (
        "initialConditions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar gammaDim
(
    "gamma",
    dimless,
    initialConditions.lookup("gamma")
);


dimensionedScalar T0Dim
(
    "T0",
    dimTemperature,
    initialConditions.lookup("T0")
);

dimensionedScalar p0Dim
(
    "p0",
    dimPressure,
    initialConditions.lookup("p0")
);


dimensionedScalar RDim
(
    "R",
    Foam::dimGasConstant,
    initialConditions.lookup("R")
);


dimensionedScalar LDim
(
    "L",
    Foam::dimLength,
    initialConditions.lookup("L")
);


dimensionedScalar reDim
(
    "re",
    Foam::dimLength,
    initialConditions.lookup("re")
);


Info << p0Dim << endl;
Info << T0Dim << endl;
Info << RDim << endl;
Info << gammaDim << endl;
Info << LDim << endl;
Info << reDim << endl;


scalar gamma=gammaDim.value();
scalar T0 = T0Dim.value();
scalar p0 = p0Dim.value(); 
scalar R = RDim.value(); 
scalar L = LDim.value(); 
scalar re = reDim.value(); 


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

Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


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


volScalarField distWall// = wallDist(mesh).y(); 
(
    IOobject
    (
        "distWall",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("distWall", dimLength, 0)
);

distWall = wallDist(mesh).y();
