
//---------------------------------------------------------------------
//
// Program to generate an initial condition file for lammps with a
// loop of DNA in a circular formation with ellipsoid atoms and a
// specified linking number.
//
//---------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#define PI 3.14159265358979

using namespace std;

struct myenums {
  // set up the types for atoms, bonds, angle
  int  BEND, TORS, TORSEND,
    DNADNA,
    DNA, DPAT, NCORE, NPAT1, NPAT2, PCORE, PPAT;
  myenums() {
    BEND=1; TORS=2; TORSEND=3;
    DNADNA=1;
    DNA=1;;
  }
} TYPE;

struct evec {
  double ei,ej,ek;
  evec cross(evec a) {
    evec c;
    c.ei=ej*a.ek - ek*a.ej;
    c.ej=ek*a.ei - ei*a.ek;
    c.ek=ei*a.ej - ej*a.ei;
    return c;
  }
  double length() {
    return sqrt(ei*ei + ej*ej + ek*ek);
  }
  void make_unit() { // makes it a unit vector
    double l;
    l=length();
    ei/=l;
    ej/=l;
    ek/=l;
  }
  //evec();
  //evec(double x,double y,double z) : ei(x), ej(y), ek(z) {}
};

struct quaternion {
  double q0,q1,q2,q3;
  double norm();
  void make_quat(evec, evec, evec);
  evec xaxis();
  evec zaxis();
  void rotate(evec,double);
  void make_unit();
  quaternion mult(quaternion);
  quaternion() {};
  quaternion(double q0,double q1,double q2,double q3) : q0(q0), q1(q1), q2(q2), q3(q3) {}
private : 
  double sign(double);
};

struct atom {
  atom(double a,double b, double c) : x(a), y(b), z(c) {};
  atom() {}
  double x,y,z,density;
  double q[4];
  int id,
    type, 
    mol,
    ellipse_flag;
  void quat(quaternion quat) {
    q[0]=quat.q0;
    q[1]=quat.q1;
    q[2]=quat.q2;
    q[3]=quat.q3;
  }
  evec zaxis() {
    return ( quaternion(q[0],q[1],q[2],q[3]) ).zaxis();
  }
  void rotate(evec a,double the) {
    // pase these parameters onto the quaternion rotate function
    quaternion myq=quaternion(q[0],q[1],q[2],q[3]);
    myq.rotate(a,the);
    quat(myq);
  }
};



struct bond {
bond(int aa, int bb, int t) : a(aa), b(bb), type(t) {};
  int a,b,
    type;
};

struct angle {
angle(int aa, int bb, int cc, int t) : a(aa), b(bb), c(cc), type(t) {};
  int a,b,c,
    type;
};



void add_DNA(atom,atom,vector<atom> &,vector<bond> &,vector<angle> &, double, double,double,double);
void add_DNA(atom,atom,vector<atom> &,vector<bond> &,vector<angle> &, double [3], double, double);



int main() {


  int N,              // DNA length
    nTw,              // Number of twists to add
    dnatype;          // atom type for DNA

  double lx,ly,lz,    // box size
    density,          // density of DNA beads
    S,                // radius of circle
    dang,             // angle between DNA beads round the circle
    d_theta,          // angle to rotate each bead by to get desired twist
    theta;            // an angle

  string fn;    // file name
  ofstream ouf;

  vector<atom> atoms;
  vector<bond> bonds;
  vector<angle> angles;

  // get some parameters
 getsize:
  cout<<"Length of DNA: "<<endl; 
  cin>>N;
  if (! N>0 ) {
    cout<<"Length must be an integer >0."<<endl; 
    goto getsize;
  }

  S=double(N)/(2.0*PI);

 lbox:
  cout<<"Enter three values for size of box:"<<endl; 
  cin>>lx>>ly>>lz;
  if (lx<1.5*S || ly<1.5*S || lz<1.5*S) {
    cout<<"Box too small to fit DNA circle."<<endl; 
    goto lbox;
  }

 ntwists:
  cout<<"Enter number of twists to add to DNA loop"<<endl;
  cin>>nTw;

 getfilename:
  cout<<"Enter output file name"<<endl;
  cin>>fn;

  // set some parameters
  dnatype = 1;
  dang = 2.0*PI/double(N);
  d_theta = nTw* 2 * PI / double(N);
  density = 6.0/PI;

  // Print a message
  cout<<endl;
  cout<<"Generating a circular DNA of length "<<N<<endl;
  cout<<"in a box of size "<<lx<<"x"<<ly<<"x"<<lz<<endl;
  cout<<"and adding "<<nTw<<" twists"<<endl;
  cout<<endl;

  // set up DNA
  {
    double x,y,z;
    atom last,lastlast;

    // 1st core
    x=S; y=0.0; z=0.0;
    atoms.push_back( atom(x,y,z) );
    atoms.back().type=dnatype;
    atoms.back().id=1;
    atoms.back().mol=1;
    atoms.back().ellipse_flag=1;
    atoms.back().density=density;
    atoms.back().q[0]=1; atoms.back().q[1]=0; atoms.back().q[2]=0; atoms.back().q[3]=0;
    last=atoms.back();
    lastlast=atom(0.0,0.0,0.0);
    lastlast.id=0;

    // rest 
    for (int i=1;i<N;i++) {
      x=S*cos(i*dang); y=S*sin(i*dang); z=0.0;
      add_DNA(last,lastlast,atoms,bonds,angles,density,x,y,z); 
      lastlast=last; 
      last=atoms.back();
    }

    // complete the loop
    bonds.push_back( bond(last.id,atoms.front().id,TYPE.DNADNA) );
    angles.push_back( angle(lastlast.id,last.id,atoms[0].id,TYPE.BEND) );
    angles.push_back( angle(last.id,atoms[0].id,atoms[1].id,TYPE.BEND) );
    angles.push_back( angle(last.id,atoms[0].id,atoms[1].id,TYPE.TORS) );

  }


  // do orientations
  { 
    evec zax,xax,yax;  // vectors
    quaternion q;      // quaternion

    yax.ei=0; yax.ej=0; yax.ek=-1;

    for (int i=0;i<N-1;i++) {
      zax.ei=atoms[i+1].x-atoms[i].x; zax.ej=atoms[i+1].y-atoms[i].y; zax.ek=atoms[i+1].z-atoms[i].z;
      zax.make_unit();
      xax=yax.cross(zax);
      q.make_quat(xax,yax,zax);
      atoms[i].quat(q);
    }
    zax.ei=atoms[0].x-atoms[N-1].x; zax.ej=atoms[0].y-atoms[N-1].y; zax.ek=atoms[0].z-atoms[N-1].z;
    zax.make_unit();
    xax=yax.cross(zax);
    q.make_quat(xax,yax,zax);
    atoms[N-1].quat(q);

  } 


  // add twists
  theta=0.0;
  for (int i=0;i<N;i++) {
    // rotate by theta about the z axis
    atoms[i].rotate(atoms[i].zaxis(),theta); 
    theta+=d_theta;         // by an angle which increases each time
  }

  // output
  ouf.open( fn.c_str() );
  ouf<<" LAMMPS data file"<<endl;
  ouf<<endl;
  ouf<<" "<<atoms.size()<<" atoms"<<endl;
  ouf<<" "<<N<<" ellipsoids"<<endl;
  ouf<<" "<<bonds.size()<<" bonds"<<endl;
  ouf<<" "<<angles.size()<<" angles"<<endl;
  ouf<<endl;
  ouf<<" 1"<<" atom types"<<endl;
  ouf<<" 1 bond types"<<endl;
  ouf<<" 2"<<" angle types"<<endl;
  ouf<<endl;
  ouf<<" "<<-0.5*lx<<" "<<0.5*lx<<" xlo xhi"<<endl;
  ouf<<" "<<-0.5*ly<<" "<<0.5*ly<<" ylo yhi"<<endl;
  ouf<<" "<<-0.5*lz<<" "<<0.5*lz<<" zlo zhi"<<endl;
  ouf<<endl;
  // +++++++++ MASSES +++++++++
  ouf<<endl<<" Masses"<<endl<<endl;
  ouf<<TYPE.DNA<<" "<<1<<endl; // mass is ignored for ellipsoids
  // +++++++++ ATOMS ++++++++++
  ouf<<endl<<" Atoms"<<endl<<endl;
  for (int i=0;i<atoms.size();i++) {
    ouf<<" "<<atoms[i].id<<" "<<atoms[i].type<<" "<<atoms[i].x<<" "<<atoms[i].y<<" "<<atoms[i].z<<" "<<atoms[i].mol<<" "<<atoms[i].ellipse_flag<<" "<<atoms[i].density<<endl;
  }
  // +++++++ ELLIPSOIDS +++++++
  ouf<<endl<<" Ellipsoids"<<endl<<endl;
  for (int i=0;i<atoms.size();i++) {
    if (atoms[i].ellipse_flag==1) {
      ouf<<" "<<atoms[i].id<<" 1 1 1 "<<atoms[i].q[0]<<" "<<atoms[i].q[1]<<" "<<atoms[i].q[2]<<" "<<atoms[i].q[3]<<endl; 
    }
  }
  // +++++++ VELOCITIES +++++++
  ouf<<endl<<" Velocities"<<endl<<endl;
  for (int i=0;i<atoms.size();i++) {
    ouf<<" "<<atoms[i].id<<" 0 0 0 0 0 0"<<endl;
  }
  // +++++++++ BONDS ++++++++++
  ouf<<endl<<" Bonds"<<endl<<endl;
  for (int i=0;i<bonds.size();i++) {
    ouf<<" "<<i+1<<" "<<bonds[i].type<<" "<<bonds[i].a<<" "<<bonds[i].b<<endl;
  }
  // +++++++++ ANGLES +++++++++
  ouf<<endl<<" Angles"<<endl<<endl;
  for (int i=0;i<angles.size();i++) {
    if (angles[i].type==TYPE.TORS && angles[i].c>atoms.size()) {angles[i].c-=atoms.size();}
    ouf<<" "<<i+1<<" "<<angles[i].type<<" "<<angles[i].a<<" "<<angles[i].b<<" "<<angles[i].c<<endl;
  }
  ouf.close();
  
  // clean up vectors
  atoms.clear(); vector<atom>().swap(atoms);
  bonds.clear(); vector<bond>().swap(bonds);
  angles.clear(); vector<angle>().swap(angles);


}


double quaternion::sign(double x) {return (x >= 0.0) ? +1.0 : -1.0;}
double quaternion::norm() { return sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3); }

void quaternion::make_quat(evec xax, evec yax, evec zax) {

  double r11,r12,r13,
    r21,r22,r23,
    r31,r32,r33,
    r;

  r11=xax.ei; r21=xax.ej; r31=xax.ek;
  r12=yax.ei; r22=yax.ej; r32=yax.ek;
  r13=zax.ei; r23=zax.ej; r33=zax.ek;

  q0=( r11 + r22 + r33 + 1.0f) / 4.0;
  q1=( r11 - r22 - r33 + 1.0f) / 4.0;
  q2=(-r11 + r22 - r33 + 1.0f) / 4.0;
  q3=(-r11 - r22 + r33 + 1.0f) / 4.0;

  if(q0 < 0.0) q0 = 0.0;
  if(q1 < 0.0) q1 = 0.0;
  if(q2 < 0.0) q2 = 0.0;
  if(q3 < 0.0) q3 = 0.0;

  q0 = sqrt(q0);
  q1 = sqrt(q1);
  q2 = sqrt(q2);
  q3 = sqrt(q3);

  if(q0 >= q1 && q0 >= q2 && q0 >= q3) {
    q0 *= +1.0f;
    q1 *= sign(r32 - r23);
    q2 *= sign(r13 - r31);
    q3 *= sign(r21 - r12);
  } else if(q1 >= q0 && q1 >= q2 && q1 >= q3) {
    q0 *= sign(r32 - r23);
    q1 *= 1.0;
    q2 *= sign(r21 + r12);
    q3 *= sign(r13 + r31);
  } else if(q2 >= q0 && q2 >= q1 && q2 >= q3) {
    q0 *= sign(r13 - r31);
    q1 *= sign(r21 + r12);
    q2 *= 1.0;
    q3 *= sign(r32 + r23);
  } else if(q3 >= q0 && q3 >= q1 && q3 >= q2) {
    q0 *= sign(r21 - r12);
    q1 *= sign(r31 + r13);
    q2 *= sign(r32 + r23);
    q3 *= 1.0;
  } else {
    cout<<"quaternion error"<<endl;;
  }
  r = norm();
  q0 /= r;
  q1 /= r;
  q2 /= r;
  q3 /= r;

}

evec quaternion::zaxis() {
  // gives the unit vector corresponding to the z-axis of the bead
  evec p;
  p.ei = 2.0*q1*q3 + 2.0*q0*q2;
  p.ej = 2.0*q2*q3 - 2.0*q0*q1;
  p.ek = q0*q0 - q1*q1 - q2*q2 + q3*q3;
  return p;
}

void quaternion::make_unit() {
  // make it a unit quaternion
  double l=norm();
  q0/=l;
  q1/=l;
  q2/=l;
  q3/=l;
}

quaternion quaternion::mult(quaternion b) {
  // multiply quaterion a by b to get c
  quaternion c;
  c.q0 = q0*b.q0 - q1*b.q1 - q2*b.q2 - q3*b.q3;
  c.q1 = q0*b.q1 + q1*b.q0 + q2*b.q3 - q3*b.q2;
  c.q2 = q0*b.q2 - q1*b.q3 + q2*b.q0 + q3*b.q1;
  c.q3 = q0*b.q3 + q1*b.q2 - q2*b.q1 + q3*b.q0;
  return c;
}


void quaternion::rotate(evec axis, double angle) {
  // rotate the quaternion by angle about axis

  quaternion b,c;
  double sinhangle;
  
  // generate a quaternion for the rotation
  b.q0 = cos( 0.5*angle );
  sinhangle = sin( 0.5*angle );
  b.q1 = sinhangle * axis.ei;
  b.q2 = sinhangle * axis.ej;
  b.q3 = sinhangle * axis.ek;

  // multiply the two quaternions
  c = b.mult(*this);

  // make sure it is still a unit quaternion
  c.make_unit();
  *this = c;

}


void add_DNA(atom last,atom lastlast,vector<atom> &atoms,vector<bond> &bonds,vector<angle> &angles, double box[3], double de, double density) {
  // Add a DNA bead to the random configuration

  double theta,phi,
    dx,dy,dz,
    hbox[3],id;

  for (int i=0;i<3;i++) {
    hbox[i]=box[i]*0.5;
  }

  // position
  do {
    theta=double(rand())/double(RAND_MAX)*PI;
    phi=double(rand())/double(RAND_MAX)*2.0*PI;
    dx=last.x+sin(theta)*cos(phi);
    dy=last.y+sin(theta)*sin(phi);
    dz=last.z+cos(theta);
  } while (abs(dx)>hbox[0]||abs(dy)>hbox[1]||abs(dz)>hbox[2]); // reject if outside box

  id=atoms.back().id;


  id++;
  atoms.push_back( atom(dx,dy,dz) );
  atoms.back().type=TYPE.DNA;
  atoms.back().id=id;
  atoms.back().mol=last.mol;
  atoms.back().ellipse_flag=1;
  atoms.back().density=density;   

  // oreintation
  atoms.back().q[0]=1; atoms.back().q[1]=0; atoms.back().q[2]=0; atoms.back().q[3]=0;

  // bond
  bonds.push_back( bond(last.id,atoms.back().id,TYPE.DNADNA) );

  // angle
  if (lastlast.id==0) { // this is bead number 2
    angles.push_back( angle(last.id,atoms.back().id,atoms.back().id,TYPE.TORS) ); // the third id here does nothing
  } else { 
    angles.push_back( angle(lastlast.id,last.id,atoms.back().id,TYPE.BEND) );
    angles.push_back( angle(last.id,atoms.back().id,atoms.back().id,TYPE.TORS) ); // the third id here does nothing   
  } 
   

}

void add_DNA(atom last,atom lastlast,vector<atom> &atoms,vector<bond> &bonds,vector<angle> &angles, double density, double x,double y,double z) {
  // Add a DNA bead to the loop configuration

  atoms.push_back( atom(x,y,z) );
  atoms.back().type=TYPE.DNA;
  atoms.back().id=last.id+1;
  atoms.back().mol=last.mol;
  atoms.back().ellipse_flag=1;
  atoms.back().density=density;  

  // oreintation
  atoms.back().q[0]=1; atoms.back().q[1]=0; atoms.back().q[2]=0; atoms.back().q[3]=0;

  // bond
  bonds.push_back( bond(last.id,atoms.back().id,TYPE.DNADNA) );

  // angle
  if (lastlast.id==0) { // this is bead number 2
    angles.push_back( angle(last.id,atoms.back().id,atoms.back().id+1,TYPE.TORS) ); // the third id here does nothing
  } else { 
    angles.push_back( angle(lastlast.id,last.id,atoms.back().id,TYPE.BEND) );
    angles.push_back( angle(last.id,atoms.back().id,atoms.back().id+1,TYPE.TORS) ); // the third id here does nothing   
  } 


}
