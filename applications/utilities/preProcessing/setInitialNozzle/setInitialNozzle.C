#include "fvCFD.H"
#include "wallDist.H"
#include "NewtonSecantRoot.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myFunction
{
private:
  scalar areaRatio_;
  scalar gamma_;
public:
  myFunction(scalar aR, scalar gamma){areaRatio_=aR; gamma_=gamma;};
//  ~myFunction();
  scalar operator()(scalar x) const { 
  return 1.0/x*(pow(2.0/(gamma_+1.0)*(1.0+((gamma_-1.0)/2.0)*pow(x,2.0)),(gamma_+1.0)/(2*(gamma_-1.0))))-areaRatio_;}
};

class myDerivative
{
private:
  scalar areaRatio_;
  scalar gamma_;
public:
  myDerivative(scalar aR, scalar gamma){areaRatio_=aR; gamma_=gamma;};
//  ~myDerivative();
  scalar operator()(scalar x) const {
/*
  scalar term1= 1.0/x;
  //scalar term3=(2.0/(gamma_+1.0))*(1.0+((gamma_-1.0)/2.0)*pow(x,2.0));
  //scalar term2= pow(term3,(gamma_+1.0)/(2.0*(gamma_-1.0)));
  scalar term2=pow((2/(gamma_+1))*(1+((gamma_-1)/2)*pow(x,2)),(gamma_+1)/(2*(gamma_-1)));
  scalar dTerm1_dM = -1.0/pow(x,2.0);
  scalar dTerm2_dM = (gamma_+1)/(2*(gamma_-1))*pow(term2,3-gamma_/(2*(gamma_-1)))*((gamma_-1)/(gamma_+1)*2*x);
  //scalar dTerm2_dM = (gamma_+1.0)/(2.0*(gamma_-1.0))*pow(term3,(3.0-gamma_)/(2.0*(gamma_-1.0)))*((gamma_-1.0)/(gamma_+1.0)*2.0*x);
  return dTerm1_dM * term2 + term1 * dTerm2_dM;
  //return dTerm1_dM * term2 + pow(term3, (3.0-gamma_)/(2.0*(gamma_-1.0)));
*/  

  scalar term1=(2.0/(gamma_+1.0))*(1.0+((gamma_-1.0)/2.0)*pow(x,2.0));
  scalar term2=pow(term1,(gamma_+1.0)/(2.0*(gamma_-1.0)));
  return term2*(1.0/term1-1.0/pow(x,2.0));

  }
};
}

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
#include "createFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//scalar rt = 35.1;
scalar eps = 1e-3;
scalar sub = 0.1;
scalar sup = 5.0;
scalar xt;
scalar rmin=GREAT;
/* // Small Test
const myFunction    f(1.2);
const myDerivative  df(1.2);

NewtonSecantRoot<myFunction,myDerivative> findRoot(f,df,1e-3);

Info << f(2) << endl;
Info << df(2) << endl;

scalar root=findRoot.root(0.5);

Info << "root = " << root << endl;
*/

forAll(mesh.C(),celli)
{
scalar x = mesh.C()[celli].x();
scalar y = mesh.C()[celli].y();
    if(x<L && y<re)
    {
       scalar r = y + distWall[celli]; 
       if(r<=rmin)
       {
          xt=x;
          rmin=r;
       }
    }
}

forAll(mesh.C(), celli)
{
scalar x = mesh.C()[celli].x();
scalar y = mesh.C()[celli].y();
scalar r = y + distWall[celli];
    if (x<L && y < re)
    {
      //scalar areaRatio = pow(r/rt,2);
      scalar areaRatio = pow(r/rmin,2);
      myFunction  f(areaRatio, gamma);
      myDerivative df(areaRatio, gamma);
      scalar M=1.0;
      
      NewtonSecantRoot<myFunction,myDerivative> findRoot(f,df,eps);
      if(x<=xt){
            M=findRoot.root(sub);
          }
      else{
            M=findRoot.root(sup);
          }
      scalar Ttemp=T0/(1+(gamma-1)/2*pow(M,2));
      scalar pTemp=p0*Foam::pow(1+(gamma-1)/2*Foam::pow(M,2),-gamma/(gamma-1));
      T[celli]= Ttemp;
      p[celli]= pTemp;
    }
}

label patchLabel = mesh.boundaryMesh().findPatchID("inlet");

Info << patchLabel << endl;
  
const fvPatch& patch = mesh.boundary()[patchLabel];
scalarField& Tsurf=T.boundaryFieldRef()[patchLabel];
scalarField& pSurf=p.boundaryFieldRef()[patchLabel];
scalarField& wallSurf=distWall.boundaryFieldRef()[patchLabel];

forAll(patch,facei)
{
scalar x = patch.Cf()[facei].x();
scalar y = patch.Cf()[facei].y();
scalar r = y + wallSurf[facei];
//scalar areaRatio = pow(r/rt,2);
scalar areaRatio = pow(r/rmin,2);
myFunction  f(areaRatio, gamma);
myDerivative df(areaRatio, gamma);
scalar M=1.0;

NewtonSecantRoot<myFunction,myDerivative> findRoot(f,df,eps);
if(x<=xt){
      M=findRoot.root(sub);
    }
else{
      M=findRoot.root(sup);
    }
scalar Ttemp=T0/(1+(gamma-1)/2*pow(M,2));
scalar pTemp=p0*Foam::pow(1+(gamma-1)/2*Foam::pow(M,2),-gamma/(gamma-1));
Tsurf[facei]= Ttemp;
pSurf[facei]= pTemp;
Info << "faceI = " << facei << " y = " << y << " r = " << r << " M = " << M << " T = " << Ttemp << " p = " << pTemp <<endl;
}

/*
patchLabel = mesh.boundaryMesh().findPatchID("Outlet");

const fvPatch& patch2 = mesh.boundary()[patchLabel];
scalarField& Tsurf2=T.boundaryFieldRef()[patchLabel];
scalarField& pSurf2=p.boundaryFieldRef()[patchLabel];
scalarField& wallSurf2=distWall.boundaryFieldRef()[patchLabel];

forAll(patch2,facei)
{
scalar x = patch2.Cf()[facei].x();
scalar y = patch2.Cf()[facei].y();
scalar r = y + wallSurf2[facei];
//scalar areaRatio = pow(r/rt,2);
scalar areaRatio = pow(r/rmin,2);
myFunction  f(areaRatio, gamma);
myDerivative df(areaRatio, gamma);
scalar M=1.0;

NewtonSecantRoot<myFunction,myDerivative> findRoot(f,df,eps);
if(x<=xt){
      M=findRoot.root(sub);
    }
else{
      M=findRoot.root(sup);
    }
//M=Foam::sqrt(M);
scalar Ttemp=T0/(1+(gamma-1)/2*pow(M,2));
scalar pTemp=p0*Foam::pow(1+(gamma-1)/2*Foam::pow(M,2),-gamma/(gamma-1));
Tsurf2[facei]= Ttemp;
pSurf2[facei]= pTemp;
Info << "faceI = " << facei << " y = " << y << " r = " << r << " M = " << M << " T = " << Ttemp << " p = " << pTemp <<endl;
}
*/

Info<< "Writing field T\n" << endl;
T.write();
Info<< "Writing field p\n" << endl;
p.write();
Info << "Writing wallDist" << endl;
distWall.write();

Info << "xt = " << xt << endl;
//Info << "rt = " << rt << endl;
Info << "rmin = " << rmin << endl;

Info << "\n end\n";
//
return(0);
}

