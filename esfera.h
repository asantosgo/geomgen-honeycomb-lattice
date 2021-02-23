#ifndef ESFERA_H
#define ESFERA_H

// Author: Arturo Santos Gómez
// Update: May 4, 2018.
// MIT License

#include <math.h> // utiliza math.h

#ifndef pi           // utiliza definición de pi
#define pi 3.141592
#endif


class Sphere {              // al parecer les gusta que las clases  empiecen con letras mayúsculas
private:                     // si no pudiera acceder a estas variables desde clases heredadas, tal vez deba ponerlas como protected
    float xPos, yPos, zPos;  // coordinates
    float D, R; 
public:
    void setPosition(float, float, float);
    void setX(float);
    void setY(float);
    void setZ(float);
    void setDiameter(float);
    // metodos get
    float X();
    float Y();
    float Z();
    float radius();
    float diameter();
    float radiusSquare();
    float volume();
};

// métodos  de asignación
void Sphere::setPosition(float x, float y, float z) { xPos=x;   yPos=y;   zPos=z; }
void Sphere::setDiameter(float diam) { D = diam; R = D/2; }
void Sphere::setX(float x) { xPos=x; }
void Sphere::setY(float y) { yPos=y; }
void Sphere::setZ(float z) { zPos=z; }
// métodos get
float Sphere::diameter(){ return D; }
float Sphere::radius(){ return R; }
float Sphere::X(){ return xPos; }
float Sphere::Y(){ return yPos; }
float Sphere::Z(){ return zPos; }
float Sphere::radiusSquare() { return pow(R,2);     }
float Sphere::volume()       { return 4*pi*pow(R,2)/3; } 
/**************************************************************/

#endif  /* ESFERA_H */