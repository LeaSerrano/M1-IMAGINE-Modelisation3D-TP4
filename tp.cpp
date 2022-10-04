// -------------------------------------------
// gMini : a minimal OpenGL/GLUT application
// for 3D graphics.
// Copyright (C) 2006-2008 Tamy Boubekeur
// All rights reserved.
// -------------------------------------------

// -------------------------------------------
// Disclaimer: this code is dirty in the
// meaning that there is no attention paid to
// proper class attribute access, memory
// management or optimisation of any kind. It
// is designed for quick-and-dirty testing
// purpose.
// -------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <GL/glut.h>
#include <float.h>
#include "src/Vec3.h"
#include "src/Camera.h"

enum DisplayMode{ WIRE=0, SOLID=1, LIGHTED_WIRE=2, LIGHTED=3 };

struct Triangle {
    inline Triangle () {
        v[0] = v[1] = v[2] = 0;
    }
    inline Triangle (const Triangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
    }
    inline Triangle (unsigned int v0, unsigned int v1, unsigned int v2) {
        v[0] = v0;   v[1] = v1;   v[2] = v2;
    }
    unsigned int & operator [] (unsigned int iv) { return v[iv]; }
    unsigned int operator [] (unsigned int iv) const { return v[iv]; }
    inline virtual ~Triangle () {}
    inline Triangle & operator = (const Triangle & t) {
        v[0] = t.v[0];   v[1] = t.v[1];   v[2] = t.v[2];
        return (*this);
    }
    // membres indices des sommets du triangle:
    unsigned int v[3];
};


struct Mesh {
    std::vector< Vec3 > vertices; //array of mesh vertices positions
    std::vector< Vec3 > normals; //array of vertices normals useful for the display
    std::vector< Triangle > triangles; //array of mesh triangles
    std::vector< Vec3 > triangle_normals; //triangle normals to display face normals

    //Compute face normals for the display
    void computeTrianglesNormals(){

        //A faire : implémenter le calcul des normales par face
        //Attention commencer la fonction par triangle_normals.clear();
        triangle_normals.clear();

        //Iterer sur les triangles
        for (unsigned int i = 0; i < triangles.size(); i++) {

            //La normale du triangle i est le resultat du produit vectoriel de deux ses arêtes e_10 et e_20 normalisé (e_10^e_20)
            //L'arete e_10 est représentée par le vecteur partant du sommet 0 (triangles[i][0]) au sommet 1 (triangles[i][1])
            //L'arete e_20 est représentée par le vecteur partant du sommet 0 (triangles[i][0]) au sommet 2 (triangles[i][2])

            Vec3 s0 = vertices[triangles[i][0]];
            Vec3 s1 = vertices[triangles[i][1]];
            Vec3 s2 = vertices[triangles[i][2]];

            Vec3 e_10 = s1 - s0;
            Vec3 e_20 = s2 - s0;

            Vec3 produitVect = Vec3::cross(e_10, e_20);

            //Normaliser et ajouter dans triangle_normales
            produitVect.normalize();

            triangle_normals.push_back(produitVect);

        }

    }

    //Compute vertices normals as the average of its incident faces normals
    void computeVerticesNormals(int weight_type){
        //Utiliser weight_type : 0 uniforme, 1 aire des triangles, 2 angle du triangle

        //A faire : implémenter le calcul des normales par sommet comme la moyenne des normales des triangles incidents
        //Attention commencer la fonction par normals.clear();
        normals.clear();

        //Initializer le vecteur normals taille vertices.size() avec Vec3(0., 0., 0.)
        normals.resize(vertices.size());

        for (unsigned int i = 0; i < vertices.size(); i++) {
            normals[i] = Vec3(0., 0., 0.);
        }

        float weight = 0;
        float P, p;
        float x_10, y_10, z_10, x_20, y_20, z_20, x_21, y_21, z_21;
        float dist_10, dist_20, dist_21;
        float a, b, c;

        //Iterer sur les triangles

            for (unsigned int i = 0; i < triangles.size(); i++) {

                //Pour chaque triangle i
                //Ajouter la normale au triangle à celle de chacun des sommets en utilisant des poids
                //0 uniforme, 1 aire du triangle, 2 angle du triangle

                if (weight_type == 0) {
                    weight = 1;
                }
               else if (weight_type == 1 || weight_type == 2) {

                    Vec3 s0 = vertices[triangles[i][0]];
                    Vec3 s1 = vertices[triangles[i][1]];
                    Vec3 s2 = vertices[triangles[i][2]];

                    /*Vec3 e_10 = s1 - s0;
                    Vec3 e_20 = s2 - s0;*/

                    x_10 = pow((s1[0]-s0[0]), 2);
                    y_10 = pow((s1[1]-s0[1]), 2);
                    z_10 = pow((s1[2]-s0[2]), 2);

                    x_21= pow((s2[0]-s1[0]), 2);
                    y_21 = pow((s2[1]-s1[1]), 2);
                    z_21 = pow((s2[2]-s1[2]), 2);

                    x_20= pow((s0[0]-s2[0]), 2);
                    y_20 = pow((s0[1]-s2[1]), 2);
                    z_20 = pow((s0[2]-s2[2]), 2);

                    dist_10 = sqrt((x_10 + y_10 + z_10));
                    dist_21 = sqrt(x_21 + y_21 + z_21);
                    dist_20 = sqrt(x_20 + y_20 + z_20);

                    if (weight_type == 1) {
                        P = dist_10 + dist_21 + dist_20;

                        p = P/2;

                        weight = sqrt(p*(p-dist_10)*(p-dist_21)*(p-dist_20));
                        // ou avec le cross product : weight = Vec3::cross(e_10, e_20)/2;
                    }

                    else {
                        a = pow(dist_10, 2);
                        b = pow(dist_21, 2);
                        c = pow(dist_20, 2);

                        if (i == triangles[i][0]) {
                            weight = acos((b + c - a) / 2*b*c);
                        }
                        else if (i == triangles[i][1]) {
                            weight = acos((a + c - b) / 2*a*c);
                        }
                        else {
                            weight = acos((a + b - c) / 2*a*b); 
                        }
                    }
                }

                for (unsigned int j = 0; j < 3; j++) {
                    float normalsX = triangle_normals[i][0]*weight;
                    float normalsY = triangle_normals[i][1]*weight;
                    float normalsZ = triangle_normals[i][2]*weight;

                    normals[triangles[i][j]] += Vec3(normalsX, normalsY, normalsZ);
                }

            }
            //Iterer sur les normales et les normaliser

            for (unsigned int k = 0; k < normals.size(); k++) {
                normals[k].normalize();
            }
    }

    void computeNormals(int weight_type){
        computeTrianglesNormals();
        computeVerticesNormals(weight_type);
    }

};

//Transformation made of a rotation and translation
struct Transformation {
    Mat3 rotation;
    Vec3 translation;
};

//Basis ( origin, i, j ,k )
struct Basis {
    inline Basis ( Vec3 const & i_origin,  Vec3 const & i_i, Vec3 const & i_j, Vec3 const & i_k) {
        origin = i_origin; i = i_i ; j = i_j ; k = i_k;
    }

    inline Basis ( ) {
        origin = Vec3(0., 0., 0.);
        i = Vec3(1., 0., 0.) ; j = Vec3(0., 1., 0.) ; k = Vec3(0., 0., 1.);
    }
    Vec3 operator [] (unsigned int ib) {
        if(ib==0) return i;
        if(ib==1) return j;
        return k;}

    Vec3 origin;
    Vec3 i;
    Vec3 j;
    Vec3 k;
};

//Fonction à completer
void collect_one_ring (std::vector<Vec3> const & i_vertices,
                       std::vector< Triangle > const & i_triangles,
                       std::vector<std::vector<unsigned int> > & o_one_ring) {//one-ring of each vertex, i.e. a list of vertices with which it shares an edge

    //Initialiser le vecteur de o_one_ring de la taille du vecteur vertices
    o_one_ring.clear();
    o_one_ring.resize(i_vertices.size());

    for (long unsigned int i = 0; i < i_vertices.size(); i++) {
        //Parcourir les triangles et ajouter les voisins dans le 1-voisinage
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                //Attention verifier que l'indice n'est pas deja present
                if (std::find(std::begin(o_one_ring[i_triangles[i][j]]), std::end(o_one_ring[i_triangles[i][j]]), i_triangles[i][k]) == std::end(o_one_ring[i_triangles[i][j]]) && i_triangles[i][j] != i_triangles[i][k]) {
                    o_one_ring[i_triangles[i][j]].push_back(i_triangles[i][k]);
                }
            }
        }

        //Tous les points opposés dans le triangle sont reliés
    }

}

//Fonction à completer
void compute_vertex_valences (const std::vector<Vec3> & i_vertices,
                              const std::vector< Triangle > & i_triangles,
                              std::vector<unsigned int> & o_valences ) {
    //Utiliser la fonction collect_one_ring pour récuperer le 1-voisinage
    std::vector<std::vector<unsigned int> > o_one_ring;
    collect_one_ring(i_vertices, i_triangles, o_one_ring);

    for (long unsigned int i = 0; i < o_one_ring.size(); i++) {
        o_valences.push_back(o_one_ring[i].size());
    }
}

// Fonctions que j'ai ajoutées
float min_valence (std::vector<unsigned int> i_valence) {
    unsigned int min = i_valence[0]; 

    for (unsigned int i = 0; i < i_valence.size(); i++) {
        if (i_valence[i]< min) {
            min = i_valence[i];
        }
    }

    return (float)min;
}

float max_valence (std::vector<unsigned int> i_valence) {
    unsigned int max = i_valence[0];

    for (unsigned int i = 0; i < i_valence.size(); i++) {
        if (i_valence[i] > max) {
            max = i_valence[i];
        }
    }

    return (float)max;
}
//


//
//TP3 Ex1

//Fonction qui permet d'afficher les courbes
void DrawCurve(std::vector<Vec3> TabPointsOfCurve, long nbPoints) {
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0f, 1.0f, 0.0f);
    for (int i = 0; i < nbPoints; i++) {
        glVertex3f(TabPointsOfCurve[i][0], TabPointsOfCurve[i][1], TabPointsOfCurve[i][2]);
    }
    glEnd();
}

//Fonction qui permet de donner une liste de points correspondants à la fonction d'Hermite
std::vector<Vec3> HermiteCubicCurve(Vec3 P0, Vec3 P1, Vec3 V0, Vec3 V1, long nbU) {
    std::vector<Vec3> P;
    P.resize(nbU);
    float F1, F2, F3, F4;
    Vec3 Pu;
    int cpt = 0;

    for (float u = 0; u < 1; u+=(1/(float)nbU)) {
        F1 = 2*pow(u, 3) - 3*pow(u, 2) + 1;
        F2 = -2*pow(u, 3) + 3*pow(u, 2);
        F3 = pow(u, 3) - 2*pow(u, 2) + u;
        F4 = pow(u, 3) - pow(u, 2);
        Pu = F1*P0 + F2*P1 + F3*V0 + F4*V1;

        P[cpt] = Pu;
        cpt++;
    }
    return P;
}

//TP3 Ex2

//Fonction factorielle
float fact(float n) {
    float factoriel = 1;
    for (float i = 1; i <= n; i++) {
        factoriel *= i;
    }
    return factoriel;
}

//Fonction qui permet de donner une liste de points correspondants à la fonction de Béziers par l'algo de Bernstein
std::vector<Vec3> BezierCurveBernstein(std::vector<Vec3> TabControlPoint, long nbControlPoint, long nbU) {
    std::vector<Vec3> P;
    P.resize(nbU);
    float Bu;
    Vec3 Pu;
    float t1, t2, t3;
    int cpt = 0;

    for (float u = 0; u < 1; u+=(1/(float)nbU)) {
        Pu = Vec3 (0., 0., 0.);
        for (int i = 0; i < nbControlPoint; i++) {
            t1 = fact(nbControlPoint-1)/(fact(i)*fact(nbControlPoint-1-i));
            t2 = pow(u, i);
            t3 = pow((1-u), (nbControlPoint-1-i));
            Bu = t1 * t2 * t3;
            Pu += Bu*TabControlPoint[i];
        }
        P[cpt] = Pu;
        cpt++;
        Pu = Vec3 (0., 0., 0.);
    }

    return P;

}

//TP3 Ex3

//Fonction qui va tracer notre courbe en prenant en paramètre une couleur
void DrawCurveByStep(std::vector<Vec3> TabPointsOfCurve, long nbPoints, Vec3 color) {
    if (nbPoints == 1) {
        glBegin(GL_POINTS);
        glColor3f(1, 0, 0);
            glVertex3f(TabPointsOfCurve[0][0], TabPointsOfCurve[0][1], TabPointsOfCurve[0][2]);
        glEnd();
    }
    else if (nbPoints > 1) {
        glBegin(GL_LINE_STRIP);
        glColor3f(color[0], color[1], color[2]);
        for (int i = 0; i < nbPoints; i++) {
            glVertex3f(TabPointsOfCurve[i][0], TabPointsOfCurve[i][1], TabPointsOfCurve[i][2]);
        }
        glEnd();
    }
}

//Fonction récursive permettant de récupérer tous les points de nos polynômes de Casteljau en fonction d'un u donné
Vec3 Casteljau (std::vector<Vec3> TabControlPoint, int r, float u, int i) {
    if (r == 0) {
        return TabControlPoint[i];
    }
    return (1-u)*Casteljau(TabControlPoint, r-1, u, i) + u*Casteljau(TabControlPoint, r-1, u, i+1);
}

//Fonction qui permet de donner une liste de points correspondants à la fonction de Béziers par les poynômes de Casteljau
std::vector<Vec3> BezierCurveByCasteljau(std::vector<Vec3> TabControlPoint, long nbControlPoint, long nbU) {
    std::vector<Vec3> P;
    P.resize(nbU);
    int cpt = 0;
    int i = 0;
    std::vector<std::vector<Vec3>> Step;
    Step.resize(4*nbU);

    for (float u = 0; u < 1; u+=(1/(float)nbU)) {
        P[cpt] = Casteljau(TabControlPoint, TabControlPoint.size()-1, u , i);
        cpt++;
    }

    //Tracé des étapes intermédiaires
    int uActuel = 1;
    int maxJ = nbControlPoint;
    int r = 0;
    Vec3 white = Vec3(1, 1, 1);

    for (int index = 0; index < 4*nbU; index++) {
        for (int j = 0; j < maxJ; j++) {
            Step[index].push_back(Casteljau(TabControlPoint, r, uActuel/(float)nbU , j));
        }
        maxJ--;
        r++;

        if (r == nbControlPoint) {
            r = 0;
            maxJ = nbControlPoint;
            uActuel++;
        }

        DrawCurveByStep(Step[index], Step[index].size(), white);
    }
   
    return P;
}

//TP4 Ex1

void DrawSurface(std::vector<Vec3> TabPointsOfSurface, long nbPointsU) {
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0f, 1.0f, 0.0f);
    for (int i = 0; i < 1; i+=(1/(float)nbPointsU)) {
        //for (int j = 0; j < 1; j+=(1/(float)nbPointsV)) {
            glVertex3f(TabPointsOfSurface[i][0], TabPointsOfSurface[i][1], TabPointsOfSurface[i][2]);
        //}
    }
    glEnd();
}

std::vector<std::vector<Vec3>> CylindricSurface (std::vector<Vec3> bezier, Vec3 P0, Vec3 P1, int nbU, int nbV) {
    std::vector<std::vector<Vec3>> P;
    P.resize(nbU);
    std::vector<Vec3> straight;
    straight.push_back(P0);
    straight.push_back(P1);
    
    Vec3 vecStraight = straight[1] - straight[0];

    Vec3 vecPas = (1./nbV) * vecStraight;

    std::vector<Vec3> test;

    DrawCurve(straight, straight.size());

    for (float v = 0; v <= nbV; v++) {
        test.clear();
        for (float u = 0; u < bezier.size(); u++) {
            test.push_back((v*vecPas) + bezier[u]);
        }
        DrawCurve(test, test.size());
    }

    std::vector<Vec3> test2;

    for (int u = 0; u < nbU; u++) {
        test2.clear();
        test2.push_back(bezier[u]);
        test2.push_back(test[u]);
        DrawCurve(test2, test2.size());
    }

    //Definir des valeurs
    return P;
}

//TP4 Ex2

//Renommer
std::vector<std::vector<Vec3>> SurfaceReglee(std::vector<Vec3> bezier1, std::vector<Vec3> bezier2, int nbU, int nbV) {
    std::vector<std::vector<Vec3>> P;
    P.resize(nbU);

    std::vector<Vec3> test;

    for (float v = 0; v < nbV; v++) {
        test.clear();
        for (float u = 0; u < nbU; u++) {
            test.push_back((1-v)*bezier1[u] + v*bezier2[u]);
        }
        DrawCurve(test, test.size());
    }

    return P;
}

//

//Input mesh loaded at the launch of the application
Mesh mesh;
std::vector< float > mesh_valence_field; //normalized valence of each vertex

Basis basis;

bool display_normals;
bool display_smooth_normals;
bool display_mesh;
bool display_basis;
DisplayMode displayMode;
int weight_type;

// -------------------------------------------
// OpenGL/GLUT application code.
// -------------------------------------------

static GLint window;
static unsigned int SCREENWIDTH = 1600;
static unsigned int SCREENHEIGHT = 900;
static Camera camera;
static bool mouseRotatePressed = false;
static bool mouseMovePressed = false;
static bool mouseZoomPressed = false;
static int lastX=0, lastY=0, lastZoom=0;
static bool fullScreen = false;

// ------------------------------------
// File I/O
// ------------------------------------
bool saveOFF( const std::string & filename ,
              std::vector< Vec3 > const & i_vertices ,
              std::vector< Vec3 > const & i_normals ,
              std::vector< Triangle > const & i_triangles,
              std::vector< Vec3 > const & i_triangle_normals ,
              bool save_normals = true ) {
    std::ofstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open()) {
        std::cout << filename << " cannot be opened" << std::endl;
        return false;
    }

    myfile << "OFF" << std::endl ;

    unsigned int n_vertices = i_vertices.size() , n_triangles = i_triangles.size();
    myfile << n_vertices << " " << n_triangles << " 0" << std::endl;

    for( unsigned int v = 0 ; v < n_vertices ; ++v ) {
        myfile << i_vertices[v][0] << " " << i_vertices[v][1] << " " << i_vertices[v][2] << " ";
        if (save_normals) myfile << i_normals[v][0] << " " << i_normals[v][1] << " " << i_normals[v][2] << std::endl;
        else myfile << std::endl;
    }
    for( unsigned int f = 0 ; f < n_triangles ; ++f ) {
        myfile << 3 << " " << i_triangles[f][0] << " " << i_triangles[f][1] << " " << i_triangles[f][2]<< " ";
        if (save_normals) myfile << i_triangle_normals[f][0] << " " << i_triangle_normals[f][1] << " " << i_triangle_normals[f][2];
        myfile << std::endl;
    }
    myfile.close();
    return true;
}

void openOFF( std::string const & filename,
              std::vector<Vec3> & o_vertices,
              std::vector<Vec3> & o_normals,
              std::vector< Triangle > & o_triangles,
              std::vector< Vec3 > & o_triangle_normals,
              bool load_normals = true )
{
    std::ifstream myfile;
    myfile.open(filename.c_str());
    if (!myfile.is_open())
    {
        std::cout << filename << " cannot be opened" << std::endl;
        return;
    }

    std::string magic_s;

    myfile >> magic_s;

    if( magic_s != "OFF" )
    {
        std::cout << magic_s << " != OFF :   We handle ONLY *.off files." << std::endl;
        myfile.close();
        exit(1);
    }

    int n_vertices , n_faces , dummy_int;
    myfile >> n_vertices >> n_faces >> dummy_int;

    o_vertices.clear();
    o_normals.clear();

    for( int v = 0 ; v < n_vertices ; ++v )
    {
        float x , y , z ;

        myfile >> x >> y >> z ;
        o_vertices.push_back( Vec3( x , y , z ) );

        if( load_normals ) {
            myfile >> x >> y >> z;
            o_normals.push_back( Vec3( x , y , z ) );
        }
    }

    o_triangles.clear();
    o_triangle_normals.clear();
    for( int f = 0 ; f < n_faces ; ++f )
    {
        int n_vertices_on_face;
        myfile >> n_vertices_on_face;

        if( n_vertices_on_face == 3 )
        {
            unsigned int _v1 , _v2 , _v3;
            myfile >> _v1 >> _v2 >> _v3;

            o_triangles.push_back(Triangle( _v1, _v2, _v3 ));

            if( load_normals ) {
                float x , y , z ;
                myfile >> x >> y >> z;
                o_triangle_normals.push_back( Vec3( x , y , z ) );
            }
        }
        else if( n_vertices_on_face == 4 )
        {
            unsigned int _v1 , _v2 , _v3 , _v4;
            myfile >> _v1 >> _v2 >> _v3 >> _v4;

            o_triangles.push_back(Triangle(_v1, _v2, _v3 ));
            o_triangles.push_back(Triangle(_v1, _v3, _v4));
            if( load_normals ) {
                float x , y , z ;
                myfile >> x >> y >> z;
                o_triangle_normals.push_back( Vec3( x , y , z ) );
            }

        }
        else
        {
            std::cout << "We handle ONLY *.off files with 3 or 4 vertices per face" << std::endl;
            myfile.close();
            exit(1);
        }
    }

}

// ------------------------------------
// Application initialization
// ------------------------------------
void initLight () {
    GLfloat light_position1[4] = {22.0f, 16.0f, 50.0f, 0.0f};
    GLfloat direction1[3] = {-52.0f,-16.0f,-50.0f};
    GLfloat color1[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat ambient[4] = {0.3f, 0.3f, 0.3f, 0.5f};

    glLightfv (GL_LIGHT1, GL_POSITION, light_position1);
    glLightfv (GL_LIGHT1, GL_SPOT_DIRECTION, direction1);
    glLightfv (GL_LIGHT1, GL_DIFFUSE, color1);
    glLightfv (GL_LIGHT1, GL_SPECULAR, color1);
    glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambient);
    glEnable (GL_LIGHT1);
    glEnable (GL_LIGHTING);
}

void init () {
    camera.resize (SCREENWIDTH, SCREENHEIGHT);
    initLight ();
    glCullFace (GL_BACK);
    glDisable (GL_CULL_FACE);
    glDepthFunc (GL_LESS);
    glEnable (GL_DEPTH_TEST);
    glClearColor (0.2f, 0.2f, 0.3f, 1.0f);
    glEnable(GL_COLOR_MATERIAL);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glMatrixMode(GL_PROJECTION);

    display_normals = false;
    display_mesh = true;
    display_smooth_normals = true;
    displayMode = LIGHTED;
    display_basis = false;
}


// ------------------------------------
// Rendering.
// ------------------------------------

void drawVector( Vec3 const & i_from, Vec3 const & i_to ) {

    glBegin(GL_LINES);
    glVertex3f( i_from[0] , i_from[1] , i_from[2] );
    glVertex3f( i_to[0] , i_to[1] , i_to[2] );
    glEnd();
}

void drawAxis( Vec3 const & i_origin, Vec3 const & i_direction ) {

    glLineWidth(4); // for example...
    drawVector(i_origin, i_origin + i_direction);
}

void drawReferenceFrame( Vec3 const & origin, Vec3 const & i, Vec3 const & j, Vec3 const & k ) {

    glDisable(GL_LIGHTING);
    glColor3f( 0.8, 0.2, 0.2 );
    drawAxis( origin, i );
    glColor3f( 0.2, 0.8, 0.2 );
    drawAxis( origin, j );
    glColor3f( 0.2, 0.2, 0.8 );
    drawAxis( origin, k );
    glEnable(GL_LIGHTING);

}

void drawReferenceFrame( Basis & i_basis ) {
    drawReferenceFrame( i_basis.origin, i_basis.i, i_basis.j, i_basis.k );
}

typedef struct {
    float r;       // ∈ [0, 1]
    float g;       // ∈ [0, 1]
    float b;       // ∈ [0, 1]
} RGB;



RGB scalarToRGB( float scalar_value ) //Scalar_value ∈ [0, 1]
{
    RGB rgb;
    float H = scalar_value*360., S = 1., V = 0.85,
            P, Q, T,
            fract;

    (H == 360.)?(H = 0.):(H /= 60.);
    fract = H - floor(H);

    P = V*(1. - S);
    Q = V*(1. - S*fract);
    T = V*(1. - S*(1. - fract));

    if      (0. <= H && H < 1.)
        rgb = (RGB){.r = V, .g = T, .b = P};
    else if (1. <= H && H < 2.)
        rgb = (RGB){.r = Q, .g = V, .b = P};
    else if (2. <= H && H < 3.)
        rgb = (RGB){.r = P, .g = V, .b = T};
    else if (3. <= H && H < 4.)
        rgb = (RGB){.r = P, .g = Q, .b = V};
    else if (4. <= H && H < 5.)
        rgb = (RGB){.r = T, .g = P, .b = V};
    else if (5. <= H && H < 6.)
        rgb = (RGB){.r = V, .g = P, .b = Q};
    else
        rgb = (RGB){.r = 0., .g = 0., .b = 0.};

    return rgb;
}

void drawSmoothTriangleMesh( Mesh const & i_mesh , bool draw_field = false ) {
    glBegin(GL_TRIANGLES);
    for(unsigned int tIt = 0 ; tIt < i_mesh.triangles.size(); ++tIt) {

        for(unsigned int i = 0 ; i < 3 ; i++) {
            const Vec3 & p = i_mesh.vertices[i_mesh.triangles[tIt][i]]; //Vertex position
            const Vec3 & n = i_mesh.normals[i_mesh.triangles[tIt][i]]; //Vertex normal

            if( draw_field && mesh_valence_field.size() > 0 ){
                RGB color = scalarToRGB( mesh_valence_field[i_mesh.triangles[tIt][i]] );
                glColor3f( color.r, color.g, color.b );
            }
            glNormal3f( n[0] , n[1] , n[2] );
            glVertex3f( p[0] , p[1] , p[2] );
        }
    }
    glEnd();

}

void drawTriangleMesh( Mesh const & i_mesh , bool draw_field = false  ) {
    glBegin(GL_TRIANGLES);
    for(unsigned int tIt = 0 ; tIt < i_mesh.triangles.size(); ++tIt) {
        const Vec3 & n = i_mesh.triangle_normals[ tIt ]; //Triangle normal
        for(unsigned int i = 0 ; i < 3 ; i++) {
            const Vec3 & p = i_mesh.vertices[i_mesh.triangles[tIt][i]]; //Vertex position

            if( draw_field ){
                RGB color = scalarToRGB( mesh_valence_field[i_mesh.triangles[tIt][i]] );
                glColor3f( color.r, color.g, color.b );
            }
            glNormal3f( n[0] , n[1] , n[2] );
            glVertex3f( p[0] , p[1] , p[2] );
        }
    }
    glEnd();

}

void drawMesh( Mesh const & i_mesh , bool draw_field = false ){
    if(display_smooth_normals)
        drawSmoothTriangleMesh(i_mesh, draw_field) ; //Smooth display with vertices normals
    else
        drawTriangleMesh(i_mesh, draw_field) ; //Display with face normals
}

void drawVectorField( std::vector<Vec3> const & i_positions, std::vector<Vec3> const & i_directions ) {
    glLineWidth(1.);
    for(unsigned int pIt = 0 ; pIt < i_directions.size() ; ++pIt) {
        Vec3 to = i_positions[pIt] + 0.02*i_directions[pIt];
        drawVector(i_positions[pIt], to);
    }
}

void drawNormals(Mesh const& i_mesh){

    if(display_smooth_normals){
        drawVectorField( i_mesh.vertices, i_mesh.normals );
    } else {
        std::vector<Vec3> triangle_baricenters;
        for ( const Triangle& triangle : i_mesh.triangles ){
            Vec3 triangle_baricenter (0.,0.,0.);
            for( unsigned int i = 0 ; i < 3 ; i++ )
                triangle_baricenter += i_mesh.vertices[triangle[i]];
            triangle_baricenter /= 3.;
            triangle_baricenters.push_back(triangle_baricenter);
        }

        drawVectorField( triangle_baricenters, i_mesh.triangle_normals );
    }
}

//Draw fonction
void draw () {



    if(displayMode == LIGHTED || displayMode == LIGHTED_WIRE){

        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_LIGHTING);

    }  else if(displayMode == WIRE){

        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        glDisable (GL_LIGHTING);

    }  else if(displayMode == SOLID ){
        glDisable (GL_LIGHTING);
        glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);

    }

    glColor3f(0.8,1,0.8);
    drawMesh(mesh, true);

    if(displayMode == SOLID || displayMode == LIGHTED_WIRE){
        glEnable (GL_POLYGON_OFFSET_LINE);
        glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
        glLineWidth (1.0f);
        glPolygonOffset (-2.0, 1.0);

        glColor3f(0.,0.,0.);
        drawMesh(mesh, false);

        glDisable (GL_POLYGON_OFFSET_LINE);
        glEnable (GL_LIGHTING);
    }



    glDisable(GL_LIGHTING);
    if(display_normals){
        glColor3f(1.,0.,0.);
        drawNormals(mesh);
    }

    if( display_basis ){
        drawReferenceFrame(basis);
    }
    glEnable(GL_LIGHTING);


}

void changeDisplayMode(){
    if(displayMode == LIGHTED)
        displayMode = LIGHTED_WIRE;
    else if(displayMode == LIGHTED_WIRE)
        displayMode = SOLID;
    else if(displayMode == SOLID)
        displayMode = WIRE;
    else
        displayMode = LIGHTED;
}



void display () {

    glLoadIdentity ();
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply ();
    draw ();

    //Hermite
    Vec3 P0 = Vec3(0, 0, 0);
    Vec3 P1 = Vec3(2, 0, 0); 
    Vec3 V0 = Vec3(1, 1, 0);
    Vec3 V1 = Vec3(1, -1, 0); 
    int nbPoints = 5;

    /*std::vector<Vec3> hermiteCubicCurve = HermiteCubicCurve(P0, P1, V0, V1, nbPoints);
    DrawCurve(hermiteCubicCurve, nbPoints);*/

    //Tracé des points utilisés
    glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(P0[0], P0[1], P0[2]);
    glVertex3f(P1[0], P1[1], P1[2]);
    glVertex3f(V0[0], V0[1], V1[2]);
    glVertex3f(V1[0], V1[1], V1[2]);
    glEnd();


    //Bernstein
    //Exemple 1 :
    Vec3 P2 = Vec3(0, 0, 0);
    Vec3 P3 = Vec3(0, 1, 0); 
    Vec3 P4 = Vec3(1, 1, 0);
    Vec3 P5 = Vec3(1, 0, 0); 
    std::vector<Vec3> TabControlPointB;
    TabControlPointB.push_back(P2);
    TabControlPointB.push_back(P3);
    TabControlPointB.push_back(P4);
    TabControlPointB.push_back(P5);

    //Exemple 2
    /*Vec3 P16 = Vec3(0, 1, 0);
    Vec3 P17 = Vec3(1, 1, 0);
    Vec3 P18 = Vec3(0.75, 0, 0);
    std::vector<Vec3> TabControlPointB;
    TabControlPointB.push_back(P16);
    TabControlPointB.push_back(P17);
    TabControlPointB.push_back(P18);*/

    long nbControlPointB = TabControlPointB.size();
    long nbUB = 10;

    /*std::vector<Vec3> bezierCurveBernstein = BezierCurveBernstein(TabControlPointB, nbControlPointB, nbUB);
    DrawCurve(bezierCurveBernstein, nbUB);*/

    //Exemple 1 tracé de points :
    glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(P2[0], P2[1], P2[2]);
    glVertex3f(P3[0], P3[1], P3[2]);
    glVertex3f(P4[0], P4[1], P4[2]);
    glVertex3f(P5[0], P5[1], P5[2]);
    glEnd();

    //Exemple 2 tracé de points :
    /*glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(P16[0], P16[1], P16[2]);
    glVertex3f(P17[0], P17[1], P17[2]);
    glVertex3f(P18[0], P18[1], P18[2]);
    glEnd();*/


    //De Casteljau
    //Exemple 1 :
    Vec3 P6 = Vec3(0, 0, 0);
    Vec3 P7 = Vec3(0, 1, 0); 
    Vec3 P8 = Vec3(1, 1, 0);
    Vec3 P9 = Vec3(1, 0, 0); 
    std::vector<Vec3> TabControlPointC;
    TabControlPointC.push_back(P6);
    TabControlPointC.push_back(P7);
    TabControlPointC.push_back(P8);
    TabControlPointC.push_back(P9);
    long nbUC = 10;

    //Exemple 2 :
    /*Vec3 P10 = Vec3(0, 0, 0);
    Vec3 P11 = Vec3(0.25, 0.75, 0); 
    Vec3 P12 = Vec3(0.4, 0.4, 0);
    Vec3 P13 = Vec3(0.5, 0, 0); 
    Vec3 P14 = Vec3(0.75, 0, 0);
    Vec3 P15 = Vec3(1, 0.5, 0);
    std::vector<Vec3> TabControlPointC;
    TabControlPointC.push_back(P10);
    TabControlPointC.push_back(P11);
    TabControlPointC.push_back(P12);
    TabControlPointC.push_back(P13);
    TabControlPointC.push_back(P14);
    TabControlPointC.push_back(P15);
    long nbUC = 30;*/

    //long nbControlPointC = TabControlPointC.size();

    /*std::vector<Vec3> bezierCurveByCasteljau = BezierCurveByCasteljau(TabControlPointC, nbControlPointC, nbUC);
    DrawCurve(bezierCurveByCasteljau, nbUC);*/

    //Exemple 1 tracé de points :
    glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(P6[0], P6[1], P6[2]);
    glVertex3f(P7[0], P7[1], P7[2]);
    glVertex3f(P8[0], P8[1], P8[2]);
    glVertex3f(P9[0], P9[1], P9[2]);
    glEnd();    

    //Exemple 2 tracé de points :
    /*glBegin(GL_POINTS);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(P10[0], P10[1], P10[2]);
    glVertex3f(P11[0], P11[1], P11[2]);
    glVertex3f(P12[0], P12[1], P12[2]);
    glVertex3f(P13[0], P13[1], P13[2]);
    glVertex3f(P14[0], P14[1], P14[2]);
    glVertex3f(P15[0], P15[1], P15[2]);
    glEnd();*/  

    std::vector<Vec3> TabControlPointCS;
    TabControlPointCS.push_back(P6);
    TabControlPointCS.push_back(P7);
    TabControlPointCS.push_back(P8);
    TabControlPointCS.push_back(P9);

    long nbUCS = 10, nbVCS = 10;
    long nbControlPointCS = TabControlPointCS.size();

    Vec3 Point0 = Vec3(0, 0, 0);
    Vec3 Point1 = Vec3(0, 0, 1);
    std::vector<Vec3> TabStraightPoint;
    TabStraightPoint.push_back(Point0);
    TabStraightPoint.push_back(Point1);

    std::vector<Vec3> bezier = BezierCurveBernstein(TabControlPointCS, nbControlPointCS, nbUCS);
    //CylindricSurface (bezier, Point0, Point1, nbUCS, nbVCS);


    Vec3 P20 = Vec3(0, 0, 0);
    Vec3 P21 = Vec3(0, 1, 0); 
    Vec3 P22 = Vec3(1, 1, 0);
    Vec3 P23 = Vec3(1, 0, 0); 
    std::vector<Vec3> TabControlPointCS1;
    TabControlPointCS1.push_back(P20);
    TabControlPointCS1.push_back(P21);
    TabControlPointCS1.push_back(P22);
    TabControlPointCS1.push_back(P23);

    Vec3 P16 = Vec3(0, 0, 1);
    Vec3 P17 = Vec3(0, 1, 1); 
    Vec3 P18 = Vec3(1, 1, 1);
    Vec3 P19 = Vec3(1, 0, 1); 
    std::vector<Vec3> TabControlPointCS2;
    TabControlPointCS2.push_back(P16);
    TabControlPointCS2.push_back(P17);
    TabControlPointCS2.push_back(P18);
    TabControlPointCS2.push_back(P19);

    long nbUCS2 = 4, nbVCS2 = 4;

    long nbControlPointCS2 = TabControlPointCS2.size();

    std::vector<Vec3> bezier1 = BezierCurveBernstein(TabControlPointCS1, nbControlPointCS2, nbUCS2);
    std::vector<Vec3> bezier2 = BezierCurveBernstein(TabControlPointCS2, nbControlPointCS2, nbUCS2);
    
    DrawCurve(bezier1, bezier1.size());
    DrawCurve(bezier2, bezier2.size());
    SurfaceReglee(bezier1, bezier2, nbUCS2, nbVCS2);
    
    glFlush ();
    glutSwapBuffers ();
}

void idle () {
    glutPostRedisplay ();
}

// ------------------------------------
// User inputs
// ------------------------------------
//Keyboard event
void key (unsigned char keyPressed, int x, int y) {
    switch (keyPressed) {
    case 'f':
        if (fullScreen == true) {
            glutReshapeWindow (SCREENWIDTH, SCREENHEIGHT);
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }
        break;


    case 'w': //Change le mode d'affichage
        changeDisplayMode();
        break;


    case 'b': //Toggle basis display
        display_basis = !display_basis;
        break;

    case 'n': //Press n key to display normals
        display_normals = !display_normals;
        break;

    case '1': //Toggle loaded mesh display
        display_mesh = !display_mesh;
        break;

    case 's': //Switches between face normals and vertices normals
        display_smooth_normals = !display_smooth_normals;
        break;

    case '+': //Changes weight type: 0 uniforme, 1 aire des triangles, 2 angle du triangle
        weight_type ++;
        if(weight_type == 3) weight_type = 0;
        mesh.computeVerticesNormals(weight_type); //recalcul des normales avec le type de poids choisi
        break;

    default:
        break;
    }
    idle ();
}

//Mouse events
void mouse (int button, int state, int x, int y) {
    if (state == GLUT_UP) {
        mouseMovePressed = false;
        mouseRotatePressed = false;
        mouseZoomPressed = false;
    } else {
        if (button == GLUT_LEFT_BUTTON) {
            camera.beginRotate (x, y);
            mouseMovePressed = false;
            mouseRotatePressed = true;
            mouseZoomPressed = false;
        } else if (button == GLUT_RIGHT_BUTTON) {
            lastX = x;
            lastY = y;
            mouseMovePressed = true;
            mouseRotatePressed = false;
            mouseZoomPressed = false;
        } else if (button == GLUT_MIDDLE_BUTTON) {
            if (mouseZoomPressed == false) {
                lastZoom = y;
                mouseMovePressed = false;
                mouseRotatePressed = false;
                mouseZoomPressed = true;
            }
        }
    }

    idle ();
}

//Mouse motion, update camera
void motion (int x, int y) {
    if (mouseRotatePressed == true) {
        camera.rotate (x, y);
    }
    else if (mouseMovePressed == true) {
        camera.move ((x-lastX)/static_cast<float>(SCREENWIDTH), (lastY-y)/static_cast<float>(SCREENHEIGHT), 0.0);
        lastX = x;
        lastY = y;
    }
    else if (mouseZoomPressed == true) {
        camera.zoom (float (y-lastZoom)/SCREENHEIGHT);
        lastZoom = y;
    }
}


void reshape(int w, int h) {
    camera.resize (w, h);
}

// ------------------------------------
// Start of graphical application
// ------------------------------------
int main (int argc, char ** argv) {
    if (argc > 2) {
        exit (EXIT_FAILURE);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (SCREENWIDTH, SCREENHEIGHT);
    window = glutCreateWindow ("TP HAI702I");

    init ();
    glutIdleFunc (idle);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutReshapeFunc (reshape);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    key ('?', 0, 0);

    glutMainLoop ();
    return EXIT_SUCCESS;
}

