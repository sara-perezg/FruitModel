#include <vve.h>
#include <cfloat>
#include <iostream>
#include <util/parms.h>
#include <util/palette.h>
#include <util/function.h>
#include <util/vector.h>
#include <math.h>
#include <cassert>
#include <vector>
#include <string>
#include <GL/glu.h>
#include <QFile>
#include <QAction>
#include <QMenu>
#include <QKeyEvent>
#include <QMouseEvent>
#include <geometry/projection.h>
#include <geometry/intersection.h>
#include <QDebug>
#include <QString>
#include <QApplication>
#include <QSignalMapper>
#include <algorithms/tissue.h>
#include <algorithms/complex.h>
#include "complex_grid.h"
#include <storage/storage.h>
#include <storage/complex.h>
#include <util/contour.h>
#include <QTextStream>
#include "solver.h"
#include <iostream> //these two, are for printing to a file.
#include <fstream>

using geometry::Point3d;
using namespace std;
using complex_factory::Point2u;
using std::cout;
using std::endl;
using std::string;
using std::pair;
using std::max;
using geometry::Point3d;
using geometry::Point4d;
using algorithms::insert;
// voy a poner 4 en lugar de 3 porque no se porque hay un chems[3], en algunas partes y no lo quiero tocar
const int nochems = 3; //Numero inicial de chems->3, Fer RECE y RECA
typedef util::Vector<nochems, double> Chemicals;
typedef solver::Solver<nochems> Xsolver;
typedef Xsolver::tag_t tag_t;


// Plasma membrane/wall
struct PM
{
  Point3d pos;
  Point3d force;
  Point3d velocity;
  int savedId;
   Chemicals chems, derivs;
  // What should I store in XML file?
  bool serialize(storage::VVEStorage& store)
  {
    bool result = store.field("Pos", pos);
    result &= store.field("Force", force);
    result &= store.field("Velocity", velocity);
    return result;
  }
};
// Cell
struct CYT
{
  Point3d pos;
  Point2u pos_grid;
  double area;
  QString type;
  QString state;
  int savedId;
  // chemicals

  Chemicals chems, derivs;

  // required by solver class
  Xsolver::CYT_Internals interns_CYT;

  bool serialize(storage::VVEStorage& store)
  {
    bool result = store.field("Pos", pos);
	result &= store.field("chems", chems);
	result &= store.field("type", type);
	result &= store.field("state", state);
    result &= store.field("Area", area);
	// result &= store.field("Distance", distance_from_source);
    return result;
  }
};

// membrane to membrane
struct PM_PM
{
  Point3d pos;
  double area;
  Chemicals chems, derivs;

  bool serialize(storage::VVEStorage& store)
  {
    bool result = store.field("Pos", pos);
    result &= store.field("Area", area);
    return result;
  }
};

// cell to cell
struct PD
{
  Point3d pos;
  double area;
  Chemicals chems, derivs;
  bool serialize(storage::VVEStorage& store)
  {
    bool result = store.field("Pos", pos);
    result &= store.field("Area", area);
    return result;
  }
};

// cell to membrane
struct CYT_PM
{
  Point3d pos;
  double area;
  Chemicals chems, derivs;
  Xsolver::CYT_PM_Internals interns_CYT_PM;

  bool serialize(storage::VVEStorage& store)
  {
    bool result = store.field("Pos", pos);
    result &= store.field("Area", area);
    return result;
  }
};

// membrane to cell
struct PM_CYT
{
  Point3d pos;
  double area;

  Chemicals chems, derivs;

  bool serialize(storage::VVEStorage& store)
  {
    bool result = store.field("Pos", pos);
    result &= store.field("Area", area);
    return result;
  }
};



class FruitModel;

// Type of the VV graph
typedef tissue::Tissue<FruitModel,CYT,PM,PM_PM,PD,CYT_PM,PM_CYT> Tissue;

EXPORT_COMPLEX_TYPES(Tissue);

// Class defining your model
// If you change the name, don't forget to change the last line of the file too
class FruitModel : public Model
{
public:

  util::Palette palette;
  Tissue T;
  string file_ful;
  string file_fert;
  int cellFactor;
  int THeight;
  int TWidth;
  double midpoint;
  double MembraneGap;
  Matrix3d growth_tensor;
  // Solver
  Xsolver solve;
  double time;
  double time_max;
  double dt;
  double CellMaxArea;
  int speed;
  int option;
  double pF, pB, dF, pD, pBR, pRE, dRE, kF1,kF2, dRA, k; // pF:fertilization signal pB:production basal dF:production signal degradation pD:diffusion of fertilization signal
  int times;
  void readParms()
  {
    util::Parms parms("view.v");

    // read the parameters here
    T.readParms(parms, "Tissue");
    T.readViewParms(parms, "View");
	parms("Tissue","CellMaxArea",CellMaxArea);
	solve.readParms(parms, "Solver");
	parms("View","MembraneGap",MembraneGap);

	parms("Main","times",times);
	parms("Main","file_ful",file_ful); //fulox35S1_ful es el gen que se reprime por fer
	parms("Main","file_fert",file_fert); //almacenamiento de la fertilizacion

	parms("Grid","THeight",THeight);
	parms("Grid","TWidth",TWidth);
  parms("Grid","midpoint",midpoint);
	parms("Grid","cellFactor",cellFactor);


	parms("Growth","tensor",growth_tensor);


  parms("Parameters","option",option);
  parms("Parameters","speed",speed);
	parms("Parameters","pF",pF);
	parms("Parameters","pB",pB);
	parms("Parameters","dF",dF);
	parms("Parameters","pD",pD);
  // ERROR?? pBR was not delcare in this scope
  parms("Parameters","pBR",pBR);
  parms("Parameters","pRE",pRE);
  parms("Parameters","dRE",dRE);
  parms("Parameters","kF1",kF1); //why if I type KF (uppercase) is orange?
  parms("Parameters","kF2",kF2);
  parms("Parameters","dRA",dRA);
  parms("Parameters","k",k);
  }

  // Here, reread the files when they are modified
  void modifiedFiles( const std::set<std::string>& filenames )
  {
    forall(const std::string& fn, filenames)
    {
      if(fn == "pal.map")
        palette.reread();
      else if(fn == "view.v")
        readParms();
    }
  }

 // bool serialize() { return true; }

  bool serialize(storage::VVEStorage& store)
  {
    bool result = store.field("T", T);
	result &= store.field("time", time);
    return result;
  }

  // constructor
  FruitModel(QObject *parent)
    : Model(parent)
    , palette("pal.map")
    , T(palette, this)
  {
    readParms();

    // Registering the configuration files
    registerFile("pal.map");
    registerFile("view.v");

	buildTissue();
	time = 0;
	dt = 0.001;
	solve.initialize();
	SetInitials();

   }

   ~FruitModel(){ //esto es un destructor, que cuando se cierra el programa actua. esto es el destrcutor de nuestro modelo.
     //some code
   }
//-----------------------------------------------------------------
  double logistic_f(double x,double k, int sig){//, double L, float k, int x0){
    int L= 1; double x0=0.5; //float k=2;
    if (sig==-1){k=-k;}
    return (L/(1+exp(-k*(x-x0))));
  }
//-----------------------------------------------------------------
  void buildTissue(){

    // Initialize your structure here
    cell c;
    std::vector<junction> junctions;

	bool creategrid;
	creategrid = complex_factory::square_grid(THeight,TWidth,T,*this, true);

// comentet stuff for index start in (1,1)
  int i=0; //i=1
  int j=TWidth-1; //j=TWidth;

  forall(const cell& c, T.C) {
    vvcomplex::FindCenter(c, T);

    c->pos_grid[0]=i;
    c->pos_grid[1]=j;
    j--;
    // if (j==0)                       j=TWidth;
    if (j==-1){c->pos_grid[0]=i; i++; j=TWidth-1;} //reset del width, y suma de fila
  }
  ofstream myfile;
 myfile.open ("C:/Users/xfactor/Desktop/Fruit_Model/output_file.csv");
 myfile << "Writing this to a file.\n";
 myfile.close();
}
//---------------------------------------------------------------
  void calcArea()
  {
    forall(const cell& c, T.C)
    {
      vvcomplex::FindCenter(c, T);
      double a = 0;
      const Point3d& cpos = c->pos;
      forall(const junction& j, T.S.neighbors(c))
      {
        const junction& jn = T.S.nextTo(c, j);
        a += geometry::triangleArea(cpos, j->pos, jn->pos);
      }
      c->area = a;
    }
  }
//--------------------------------------------------------------------
 void SetInitials(){
	//  Initial conditions
  midpoint+=floor((THeight-1)/2);
  pF+=0.25*THeight+0.25; //esto son rectas par no tener que modificar el pD y el pF en funcion de la tHeight
  pD+=0.125*THeight-0.375;
  forall(const cell &c, T.C){
  	for (int i=0; i < nochems; i++) c->chems[i] = 1e-3;
    if (c->pos_grid[0]==midpoint){c->type= "Source";} //como el source smp va a estar en la mitad lo coloco desde ya, lo puedo  quitar al ppio asi que da un poco igual, pero asi es mas facil.
  }
}
// ------------------------------------------------------------------
  void updateDerivatives_CYT(const cell& c, const tag_t&)
  {

  double F = 0.0;
  double RECE = 0;
  double RECA = 0;

//Fertilization signal
  if (c->type == "Source"){
  	F =  pF;
  }
  else if (c->type == "Sink" ){}
//c->chems[0]=F, cogemos la primera posicion de chem, que es F.
// aqui actualizamos la ferilizacion, es la ode de fer.
// vamos anadiendo la produccion basal y quitando la degradacion
   F +=  pB - dF * c->chems[0]; //dF como porcentaje 0.1 es 10%
   if (c->chems[0]<0){c->chems[0]=0;}

  forall(const cell &d, T.C.neighbors(c)){
  //aqui estamos difundiendo la fertilizacion por el grid.
  //d es la celula vecina y c la celula origen???
  	F += pD * (d->chems[0] - c->chems[0]);
  }

  double N_pos_grid;
  int sig;

// LEFT PART OF THE GRID k<0.
  if (c->pos_grid[0]<midpoint){
    N_pos_grid=(c->pos_grid[0])/(midpoint-1);
    sig=-1;
  }
// RIGHT PART OF THE GRID k>0.
  else if (c->pos_grid[0]>midpoint) {
    N_pos_grid=(c->pos_grid[0]-(midpoint+1))/((THeight-1)-(midpoint+1));
    sig=1;
  }
  else{N_pos_grid=0; sig=1;} //el midpoint smp va a ser 0., el signo es 1. Podria ser al reves, ser N_pos_grid=1 pero sig=-1.

// hay que poner chems[], se van actualizando automaticamente por el solver

  if (option==1){
    RECE = pBR + pRE*logistic_f(N_pos_grid,k,sig) - dRE;
    RECA = (c->chems[0] / (kF1 + c->chems[0])) - dRA*c->chems[2];
  }
  else {
      RECE = pBR + pRE - dRE;
      RECA = (c->chems[0] / (kF2 + c->chems[0])) - dRA*c->chems[2];
  }


  // std::cout << c->pos_grid c->chems[1] << '\n';
   c->derivs[0]=F; // Fertilization signal
   c->derivs[1]=RECE; // Receptor Expression
   c->derivs[2]=RECA; // Receptor Activity

  }

  //------------------------------------------------------------------------------
  void updateDerivatives_CYT_PM(const cell& c, const junction& n, Tissue& T,const tag_t& )
  {
     T.S.edge(c,n)->derivs[0]=0;
     T.S.edge(c,n)->derivs[1]=0;
     T.S.edge(c,n)->derivs[2]=0;
  }

   Chemicals& derivatives_CYT(const cell& c, const tag_t& )
   {
      return c->derivs;
   }

   Chemicals& derivatives_CYT_PM(const cell& c, junction n, Tissue& T,const tag_t& )
   {
      return T.S.edge(c,n)->derivs;
   }


  Chemicals& values_CYT(const cell& c, const tag_t& )
  {
      return c->chems;
  }
  Chemicals& values_CYT_PM(const cell& c, junction n, Tissue& T,const tag_t& )
  {
      return T.S.edge(c,n)->chems;
  }


   // Methods required for SOlver class
  Xsolver::CYT_Internals& Internals_CYT(const cell& c, const tag_t& )
  {
      return c->interns_CYT;
  }

   // Methods required for SOlver class
  Xsolver::CYT_PM_Internals& Internals_CYT_PM(const cell& c, const junction& n, Tissue& T, const tag_t&  )
  {
      return T.S.edge(c,n)->interns_CYT_PM;
  }

  //----------------------------------------------------------
  void chem_step()
  {
  	// RUn Solver
  	solve(T, *this);
  }
//----------------------------------------------------------
  void step()
  {
  for (int i = 0; i < speed; i++) {
    chem_step();
  }
  /*Chemical steps*/
  // chem_step();
  }
//-----------------------------------------------------------
  void preDraw()
  {
    glClearColor(0,0,0,1);
    T.preDraw();
  }

  void postDraw()
  {
    T.postDraw();
  }
//------------------------------------------------------------
  void draw(Viewer* viewer)
  // es una funcion que se llama todo el rato, aqui se llaman a las funciones
  // drawCustomCell->citoplasma drawCellMembrane->membrana plasmatica
  // drawCircleInside-> dibuja circulos dentro de las celulas (van a ser mis RECE)
  {
    Vec vmin(HUGE_VAL, HUGE_VAL, HUGE_VAL), vmax(-HUGE_VAL,-HUGE_VAL,-HUGE_VAL);
    glLineWidth(2.0);
    glEnable(GL_POLYGON_OFFSET_FILL);
	glEnable(GL_POLYGON_OFFSET_LINE);

	forall(const junction& j, T.W)
    {
      const Point3d& jpos = j->pos;
      if(jpos.x() < vmin.x)
        vmin.x = jpos.x();
      if(jpos.y() < vmin.y)
        vmin.y = jpos.y();
      if(jpos.z() < vmin.z)
        vmin.z = jpos.z();
      if(jpos.x() > vmax.x)
        vmax.x = jpos.x();
      if(jpos.y() > vmax.y)
        vmax.y = jpos.y();
      if(jpos.z() > vmax.z)
        vmax.z = jpos.z();
    }

    forall(const cell& c, T.C)
    {
     forall(const junction& j, T.S.neighbors(c))	{
	    drawCellMembrane(c,j);
     }
	 drawCircleInside(c);
	 drawCustomCell(c);
   }
    glColor3f(1,1,1);
    viewer->setSceneBoundingBox(vmin, vmax);

    viewer->drawText(20, 20,  "Time:\t" + QString::number(time));
    // viewer->drawText(20, 20,  "Steps: " + QString::number(c->chems[]));

  }
//-----------------------------------------------------------------
  void drawCircleInside(const cell& c)
  {

  	junction fv = T.S.anyIn(c);
    junction k = fv;
    	// Draw inside circle
	   glPolygonOffset(-2, 1.0);

    // Wall endpoints
    junction li = T.S.nextTo(c, k);

	glNormal3dv(normal(c).c_data());

  // palette.blend(inicial color, final color, intensidad) si intensidad=1 vamos a tomar el color final.
  // es como un gradiente, quier que vayas de incial a final con esta intensidad.

    double CReca=c->chems[2]/10; //RECA
    if (CReca > 1) {CReca=1;}
    else if (CReca <0){CReca=0;} //si se va la concentracion a negativo sigue dibujando los circulos, me pasa con la expresion de RECA, que al tener dRA si no hay fertilizacion se va a menos infinito RECA y se dibujan circulos.
    palette.blend(1, 6, CReca); //pongo 1, en el caso de que no haya actividad en el caso de que no llegue la fertilizacion, voy a ser capaz de ver el cirlulo.
    // si hubiese receptor pero no actividad, entonces no veriamos nada y no estariamos visualizando bien el output.

    // palette.blend(0, 6, 1);
    glBegin(GL_TRIANGLE_FAN); // glBegin(lo que quiero que me dibujes) en este caso triangulos, tmb puede ser rectangulos, hexagonos etc
    glVertex3dv(c->pos.c_data());

double CRece=(c->chems[1])/(400); //RECE

// de i= numero de triangulos que quiero que me dibuje
	for (int i=0; i < 50; i++) {

      Point3d spotRe=c->pos;

  	if (CRece > 1) CRece=1;
    else if (CRece <0){CRece=0;}

    // aqui calculamos el size de los circulos, al multiplicar por la concentracion
    // lo que hacemos es aumentar o el eje x o el eje y del triangulo
    // de ahi lo de multiplicar en spot.x por el cos y en .y por el sin.
  	spotRe.x()+=(c->pos-k->pos).norm()* CRece *cos(2*M_PI*i/45)/2;
  	spotRe.y()+=(c->pos-k->pos).norm()* CRece *sin(2*M_PI*i/45)/2;

      glVertex3dv(spotRe.c_data());
  }
    glEnd();
  }
//-----------------------------------------------------------------
  void drawCustomCell(const cell& c) {
// esto es el citoplasma, que va a ser mi fertilizacion
	  // Set offset for drawing PIN
   glPolygonOffset(0.0, 1.0);

	junction fv = T.S.anyIn(c);
    junction k = fv;

    // Wall endpoints
    junction l = T.S.nextTo(c, k);

  glNormal3dv(normal(c).c_data());
  //std::cout << c->chems[0] << std::endl;

  double concentration=c->chems[0]/10; //pongo la fertilization
  // double concentration=0;
  if (concentration > 1) concentration=1;

    palette.blend(0, 1, concentration);

    glBegin(GL_TRIANGLE_FAN);

	glVertex3dv(c->pos.c_data());
    glVertex3dv(k->pos.c_data());
	glVertex3dv(l->pos.c_data());

    while(l != fv) {
      l = T.S.nextTo(c, l);
      glVertex3dv(l->pos.c_data());
    }
    glEnd();
  }
//-----------------------------------------------------------------
  void drawCellMembrane(const cell& c, const junction& k) {

	glPolygonOffset(-10.0, 1.0);
	  // Get vectors towards neighbors and normalize
      const junction& j = T.S.prevTo(c, k);
      const junction& l = T.S.nextTo(c, k);
	  const junction& ln = T.S.nextTo(c, l);


	  // corners
	  Point3d jkn=k->pos - j->pos;
	  Point3d lkn=l->pos-k->pos;
	  Point3d lnln=ln->pos-l->pos;

	  // quads
	 Point3d quadk, quadl;

	   // Find a vector to make inside points of quads
      Point3d kprev = (normal(k) ^ jkn).normalize();
      Point3d knext = (normal(k) ^ lkn).normalize();
	  Point3d lprev = (normal(l) ^ lkn).normalize();
      Point3d lnext = (normal(l) ^ lnln).normalize();


	  quadk = (kprev + knext).normalize() / sin((M_PI - acos(kprev * knext))*0.1);
      quadl = (lprev + lnext).normalize() / sin((M_PI - acos(lprev * lnext))*0.1);


	  palette.blend(0, 7, 1);

	 junction nj=T.S.nextTo(c,k);
     cell d = T.S.prevTo(k,c);
	 junction np= T.S.prevTo(d,k);

	  if (d && d != c && nj == np)
	  {

	    double concentration=c->chems[1];
		//std::cout << c->chems[1] << std::endl;

		if (concentration > 1) concentration=1;

		if (c->type == "Source" )
			palette.blend(0, 6, 1);

		//std::cout << "ALL RIGHT" << std::endl;
	  }

	 // Draw membranes

          glBegin(GL_QUADS);
          glVertex3dv(l->pos.c_data());
          glVertex3dv(k->pos.c_data());
          glVertex3dv((k->pos + MembraneGap*quadk).c_data());
          glVertex3dv((l->pos + MembraneGap*quadl).c_data());
          glEnd();
  }
//-----------------------------------------------------------------
  void drawCellWall(const cell& c, const junction& k) {
// no se la diferencia entre cell membrane y cell wall, en mi modelo,
// esto sera cuando tiene tina membrana y el cell wall distinto rollo
// pared celular hexagonal y la membrana de menor tamanio.

	  const cell& cp=T.S.prevTo(k,c);
	  const junction& l=T.S.nextTo(c,k);
      const junction& ln=T.S.nextTo(cp,k);

	  glPolygonOffset(5, 1.0);
	  palette.blend(0, 7, 1);
	 // Draw membranes

          glBegin(GL_QUADS);
          glVertex3dv(T.S.edge(c,l)->pos.c_data());
          glVertex3dv(T.S.edge(c,k)->pos.c_data());
          glVertex3dv(k->pos.c_data());
          glVertex3dv(l->pos.c_data());
          glEnd();
  }
//-----------------------------------------------------------------
  void drawWithNames()
  {
    //dibuja en la pantalla
    int i = 0;
    forall_named(const cell& c, T.S, cells)
    {
      glPushName(i++);
      glBegin(GL_POLYGON);
      forall(const junction& j, T.S.neighbors(c))
      {
        glVertex3dv(j->pos.c_data());
      }
      glEnd();
      glPopName();
    }
  }
//-----------------------------------------------------------------
  // Methods needed by the new tissue
  Point3d position(const cell& c) const { return c->pos; }
  void setPosition(const cell& c, const Point3d& p) { c->pos = p; }

  Point3d position(const junction& j) const { return j->pos; }
  void setPosition(const junction& j, const Point3d& p) { j->pos = p; }
  void setPositionHint(const junction&, const junction&, const junction&, double) {}

  Point3d normal(const cell&) const { return Point3d(0,0,1); }
  Point3d normal(const junction&) const { return Point3d(0,0,1); }

 // Grid positions

  void saveGridCellPosition(const cell& c, Point2u pos, Tissue&) {c->pos_grid = pos; }
  Point2u gridCellPosition(const cell& c, Tissue&) { return c->pos_grid;}

//-----------------------------------------------------------------
    void updateFromOld(const cell& c1, const cell& c2, const cell& c,const Tissue::division_data&, Tissue& T)
	{
		c1->chems[0]=c2->chems[0]=0.5*c->chems[0];
		c1->chems[1]=c2->chems[1]=0.5*c->chems[1];
		c1->chems[2]=c2->chems[2]=0.5*c->chems[2];
		c1->type=c2->type=c->type;
		c1->state=c2->state=c->state;

	}

};

#include <QMouseEvent>

class FruitViewer : public Viewer
{
  Q_OBJECT;

    QList<QString> cell_behaviour;
	QList<QString> cell_type;
	QList<QString> cell_state;

	QList<QAction *> cb_list;

	QList<QAction *> ct_list;

	QList<QAction *> cs_list;

	QAction *monitor_cells;
	QAction *set_initials;
	QSignalMapper* signalMapper_behave;
	QSignalMapper* signalMapper_state;

public:
  FruitViewer (QWidget* parent,const QGLWidget* shareWidget, int idx = -1)
    : Viewer(parent,shareWidget,idx)
  {

	cell_behaviour << "None" << "Source" << "WeakSource" << "StrongSource" << "Sink" << "WeakSink" << "StrongSink" ;
    cell_type << "None" << "Epidermis" << "Suspensor" << "Cortex" << "Endodermis" << "Vascular" << "Pericycle" ;
	cell_state << "None" << "Meristem" << "Elongation" << "Differentation" << "Internal" << "Surface" << "Dead" ;

	for (int i = 0; i < cell_behaviour.length(); i++) cb_list.append(new QAction(cell_behaviour.at(i), this));
	for (int i = 0; i < cell_type.length(); i++) ct_list.append(new QAction(cell_type.at(i), this));
	for (int i = 0; i < cell_state.length(); i++) cs_list.append(new QAction(cell_state.at(i), this));
	monitor_cells = new QAction("Monitor Cells",this);
	set_initials = new QAction("Initial states",this);

	signalMapper_behave = new QSignalMapper(this);
	signalMapper_state = new QSignalMapper(this);
	for (int i = 0; i < cb_list.length(); i++)
	{
		connect(cb_list[i], SIGNAL(triggered()),signalMapper_behave, SLOT(map()));
		signalMapper_behave->setMapping(cb_list[i], cell_behaviour[i]);

	}

	for (int i = 0; i < cs_list.length(); i++)
	{
		connect(cs_list[i], SIGNAL(triggered()),signalMapper_state, SLOT(map()));
		signalMapper_state->setMapping(cs_list[i], cell_state[i]);

	}

	connect (signalMapper_behave, SIGNAL(mapped(QString)), this, SLOT(ActivatedBehave(QString))) ;
	connect (signalMapper_state, SIGNAL(mapped(QString)), this, SLOT(ActivatedState(QString))) ;

  }

  ~FruitViewer() {}

public slots:

// Behaviour
  void ActivatedBehave(const QString &text)
  {
	 FruitModel* m = dynamic_cast<FruitModel*>(model());
	 int s = selectedName();

	 if(s != -1)
	 {
		const cell& c = m->T.C[s];

		int index = cell_behaviour.indexOf(text);

		c->type =  cell_behaviour.at(index);

		std::cout << (c->type).toStdString() << std::endl;

    std::cout << "-------------------\n"<< c->pos_grid[0]<<"\t" <<c->pos_grid[1] << std::endl;
    std::cout <<"RECE\t" << c->chems[1] << std::endl;
    std::cout << "FERT\t" <<c->chems[0] << std::endl;
    std::cout << "RECA\t" <<c->chems[2] << "\n-------------------\n"<<std::endl;
   }

  }
  // State


  // Behaviour
  void ActivatedState(const QString &text)
  {
	 FruitModel* m = dynamic_cast<FruitModel*>(model());
	 int s = selectedName();

	 if(s != -1)
	 {
		const cell& c = m->T.C[s];

		int index = cell_state.indexOf(text);

		c->state =  cell_state.at(index);
		std::cout << (c->state).toStdString() << std::endl;
	 }

  }


protected:

  void mousePressEvent(QMouseEvent* event)
  {
    FruitModel *m = dynamic_cast<FruitModel*>(model());


    if(event->button() == Qt::RightButton && event->modifiers() == Qt::NoButton)
    {
      updateGL();
      update();
      event->accept();
	  std::cout << "cosik" << std::endl;

	  QMenu menu( this );
	  QMenu *cb_menu;
	  QMenu *cs_menu;

	  cb_menu = menu.addMenu(tr("Cell behaviour"));
	  cs_menu = menu.addMenu(tr("Cell state"));

	  QList<QAction *>::const_iterator i;
	  for (i = cb_list.constBegin(); i != cb_list.constEnd(); ++i) cb_menu->addAction(*i);

	  menu.addSeparator();

	   for (i = cs_list.constBegin(); i != cs_list.constEnd(); ++i) cs_menu->addAction(*i);

	  menu.exec(event->globalPos());
    }
    else
    {
      Viewer::mousePressEvent(event);
    }
  }

};

#include "model.moc"


// Define 'MyModel' as the model to be used
DEFINE_MODEL(FruitModel);
DEFINE_VIEWER(FruitViewer);
