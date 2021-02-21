/**********************************************************
            GEOMETRY GENERATOR FOR DDSCAT
                = Monolayer  NPs =
(basado en el código de la monocapa de esferas dieléctricas)
   
        Description: Genera el archivo shape.dat 
                    para su uso en ddscat 

            Output file: shape.dat 

            Author: Arturo Santos Gómez

                Date: Jul 19, 2018.
              Update: Sep 8, 2018 
**********************************************************/

#include <cstdio>
#include <cmath>
#include <vector> // por si deseo usar vector en vez de new y delete[]
#include "esfera.h"  // utiliza math.h
#include <algorithm>

#ifndef pi           		  // Se recomienda usar mayusculas
#define pi 3.141592           // Definí así la constante pi porque no se como se inicializa una constante en una clase
#endif                        // aparte de que me permite no volverla a redefinir con #ifndef #endif

//const float pi = 3.141592;  // los programadores recomiendan esta forma de declarar constantes

const long SIZE = 10000000;  // Tamanio del arreglo de enteros para almacenar los dipolos
                             // Me conviene usar este tamaño porque DDSCAT maneja hasta 1 millón de dipolos

/**************************************************************/
class SphereDip : public Sphere {  // hereda las propiedades y métodos de la clase esfera
private:
    // int Ndip=0; //se puede inicializar así pero no es preferido por los programadores
    int Ndip;
   // int Rdip;
public:
   // SphereDip(float d) { Ndip=0; Rdip=R/d; }  // constructor, inicialiaza a cero 
    SphereDip() { Ndip=0; }  // constructor, inicialiaza a cero 
    void add_dipole()  { Ndip = Ndip+1; }
    int Ndipoles()     { return Ndip;   } 
  //  float volumeDip()  { return 4*pi; }
};
/**************************************************************/

int main(void){

    /************************ Declaracion de variables ************************/
    //SphereDip Sphere1, Sphere2, Sphere3, Sphere4, Sphere5;   // se crean los objetos esfera, las posiciones y los diámetros se dan mas adelante
    SphereDip Nanoparticle1, Nanoparticle2, Nanoparticle3, Nanoparticle4, Nanoparticle5, Nanoparticle6;  // se crean los objetos para las NPs

    // for geometry
    float D, R, Rbis, Vsph;              // diameter, radius, virtual radius, volume
    float dNP, r, rbis, vNP;             // for geometry of the NPs
    float H, sep, d;                     // height ot the NPs, separation between spheres, distance between dipoles
    float Vtot, aeff, aeff_dip;          
    float Py, Pz;                        // para definir la periodicidad
 

    float f1, f2;                        // para resringir la geometría
  
    float sepCC, sepSS;                  // separación entre NPs: centro-centro y superficie superficie

    /* geometría con flotantes*/
    float Rdip, Rdipsquare, Vsphdip, rdip, rdipSquare, vNPdip, Hdip;
    //float SPH1, SPH2, SPH3, SPH4, SPH5; // para definir los puntos interiores de cada esfera
    float NP1, NP2, NP3, NP4, NP5, NP6;           // para definir los puntos interiores de cada NP
    
    float X0, Y0, Z0;                   // para definir el origen de la geometría

    // coordenadas para las esferas y material del objeto
    std::vector<long> JX_data(SIZE), JY_data(SIZE), JZ_data(SIZE), ICOMP(SIZE);   // To store data, se utiliza la libreria vector. 
                                                                                  // Se usan igual que los arreglos!
                                                                                  // Es mejor usar este (que new y delete[]) 
                                                                                  // porque se pueden adaptar facilmente a vectores con datos flotantes, dobles o long

    long N, nNP, Ntot;    // for number of dipoles of one sphere: expected: N

    /* counters for number of dipoles: J stands for the total dipoles calculated */
    long J = 0;   // es importante que el tipo de dato se defina de acuerdo a SIZE

    /* iteradores del grid con enteros*/
    long JX, JY, JZ;
    long JXMIN, JYMIN, JZMIN;
    long JXMAX, JYMAX, JZMAX;
    long XDmin, XDmax, YDmin, YDmax, ZDmin, ZDmax; 

    FILE *pFile;
    pFile = fopen ("shape.dat","w");  // output file
    /**********************************************************************************/

    /********************** PARAMETROS A MODIFICAR **********************/
     D = 129.204;  /* Diameter of the sphere (in nm) */      
     dNP = 30.0;  /* radio de las NPs */                    
     sep = 0.7;   /* separación entre esferas dielectricas. Es introducida para evitar traslape entre esferas */             
     d = .7;      /* distance between dipoles (the same units than D) */
    /*********************************************************************/


  /******************************************************************************/
   R = D/2;                             // radio esfera dieléctrica
   r = dNP/2;                           // radio NP metálica
   Rbis = R + sep/2;                    // radio ficticio de la esfera dieléctrica
   rbis = r + sep/2;                    // virtual radius of the NPs
   Py = 2*Rbis;  Pz = 2*sqrt(3)*Rbis;   /* definición de la periodicidad en unidades físicas  */
   X0 = 0;  Y0 = 0;  Z0 = 0;            // Para definir el centro en shape.dat (en unidades físicas) 
  /******************************************************************************/

    sepCC=2*Rbis/sqrt(3); // separación centro-centro
    sepSS=sepCC-2*r;      // separación superficie-superficie

   /******** Cantidades que comparten las esferas ********/
   Rdip = R/d;                   // radio en términos de dipolos
   Rdipsquare = pow(Rdip, 2);    // radio al cuadrado
   rdip = r/d;
   rdipSquare = pow(rdip, 2);
   Vsph = 4*pi*pow(R,3)/3;       // volumen esfera dieléctrica
   vNP  = 4*pi*pow(r,3)/3;       // volumen NP
   N    = Vsph/pow(d,3);            // Approximation of the number of dipoles for one sphere
   nNP  = vNP/pow(d,3);            // Approximation of the number of dipoles for one sphere
  
   Ntot = 4*nNP;
   aeff = cbrt( 4*pow(r,3) );  // expected effective radius for 2 dielectric spheres and 4 NPs
  /*******************************************************/
  
   /******** Definición de la altura de las Nps ********/
   // esto no se usará, porque como no están soportadas en "nada" la altura se hace cero
   H = sqrt( pow(Rbis+rbis, 2) - pow(2*sqrt(3)*Rbis/3, 2) ); // En unidades físicas
   
   H = 0;  // por lo menciado anteriormente
   
   Hdip = H/d; // -> Hdip=0
   /****************************************************/
   
   /* Las dimensiones de las esferas se darán en términos de dipolos: */
   // en esta parte las coordenadas etan guardadas como flotantes
	// Sphere1.setX(0);   Sphere1.setY(-Rbis/d);    Sphere1.setZ(sqrt(3)*Rbis/d);    Sphere1.setDiameter(D);    // sphere 1
	// Sphere2.setX(0);   Sphere2.setY( Rbis/d);    Sphere2.setZ(sqrt(3)*Rbis/d);    Sphere2.setDiameter(D);    // sphere 2
	// Sphere3.setX(0);   Sphere3.setY(0);          Sphere3.setZ(0);                 Sphere3.setDiameter(D);    // sphere 3 (sphere in the origin)
	// Sphere4.setX(0);   Sphere4.setY(-Rbis/d);    Sphere4.setZ(-sqrt(3)*Rbis/d);   Sphere4.setDiameter(D);    // sphere 4
	// Sphere5.setX(0);   Sphere5.setY( Rbis/d);    Sphere5.setZ(-sqrt(3)*Rbis/d);   Sphere5.setDiameter(D);    // sphere 5

    // Metallic nanoparticles
    // Hdip -> 0 porque las esferas dieléctricas son virtuales.
    Nanoparticle1.setX(-Hdip);   Nanoparticle1.setY(0);          Nanoparticle1.setZ( 2*sqrt(3)*Rbis/(3*d));    Nanoparticle1.setDiameter(dNP);    // NP 1, la de arriba

	 Nanoparticle2.setX(-Hdip);   Nanoparticle2.setY(-Rbis/d);    Nanoparticle2.setZ( sqrt(3)*Rbis/(3*d));      Nanoparticle2.setDiameter(dNP);    // NP 2, izquierda arriba
	 Nanoparticle3.setX(-Hdip);   Nanoparticle3.setY(-Rbis/d);    Nanoparticle3.setZ(-sqrt(3)*Rbis/(3*d));      Nanoparticle3.setDiameter(dNP);    // NP 3, izquierda abajo

    Nanoparticle4.setX(-Hdip);   Nanoparticle4.setY(0);          Nanoparticle4.setZ(-2*sqrt(3)*Rbis/(3*d));    Nanoparticle4.setDiameter(dNP);    // NP 4, la de abajo

    Nanoparticle5.setX(-Hdip);   Nanoparticle5.setY( Rbis/d);    Nanoparticle5.setZ(-sqrt(3)*Rbis/(3*d));      Nanoparticle5.setDiameter(dNP);    // NP 5, derecha abajo 
    Nanoparticle6.setX(-Hdip);   Nanoparticle6.setY( Rbis/d);    Nanoparticle6.setZ( sqrt(3)*Rbis/(3*d));      Nanoparticle6.setDiameter(dNP);    // NP 6, derecha arriba

   printf("===================================================================\n");
   printf("\t MONOCAPA DE NPS \n");
   printf("===================================================================\n");
   printf("\t DATOS INTRODUCIDOS\n");
   printf("-------------------------------------------------------------------\n");
   printf("Diameter of the virtual dielectric sphere:    Dvirt = %.3f nm\n", D);
   printf("Diameter of the NP:                             dNP = %.3f nm\n", dNP);
   printf("Separación entre esferas virtuales:             sep = %.3f nm\n", sep);
   printf("Separación entre NPs (centro-centro):         sepCC = %.3f = %.2fr \n", sepCC, sepCC/r);
   printf("Separación entre NPs (superficie-superficie): sepSS = %.3f = %.2fr \n", sepSS, sepSS/r);
   printf("Distance between dipoles:                         d = %.3f nm\n\n", d);
   printf("Radius of the sphere:                             R = %.3f nm\n", R);
   printf("Radio de la NP:                                   r = %.3f nm\n", r);
   printf("La altura de las NPs es:                          H = %.3f nm\n", H);
   printf("Centro de masa:                                 (X0, Y0, Z0) = (%.2f nm, %.2f nm, %.2f nm)  \n", X0*d, Y0*d, Z0*d);
   printf("Centro de masa (en términos de dipolos):        (X0, Y0, Z0) = (%.2f, %.2f, %.2f)  \n", X0, Y0, Z0);
   printf("===================================================================\n");
   printf("\t EXPECTED VALUES\n");
   printf("-------------------------------------------------------------------\n");
   printf("Radius of the sphere in terms of dipoles:    Rdip = %.2f\n", Rdip);
   printf("Radius of the NP in terms of dipoles:        rdip = %.2f\n", rdip);
   printf("Altura de las NPs en dipolos:                Hdip = %.3f\n", Hdip);
   printf("Dipoles for the sphere:                        N' = %d\n", N);
   printf("Dipoles for the NP:                          nNP' = %d\n", nNP);
   printf("Volumen de la esfera:                       Vsph' = %1.3E nm^3\n", Vsph);
   printf("Volume of the NP:                            vNP' = %1.3E nm^3\n\n", vNP);

   printf("Número total de dipolos                     Ntot' = %d\n", Ntot);
   printf("Efective radius:                            aeff' = %.2f nm = %.5f um\n", aeff, aeff/1000);


  /******** definiendo rango para el mallado **********/
  JXMIN = -2*Rbis/d;    JXMAX = 2*Rbis/d;
  //JXMIN = 0;    JXMAX = D/d;
  // JXMIN = -D/d;    JXMAX = 0;
  JYMIN = -2*Rbis/d;    JYMAX = 2*Rbis/d;
  JZMIN = -3*Rbis/d;    JZMAX = 3*Rbis/d;
  /*****************************************************/

 /*************************************************************************************/
  /*************************************************************************************/
 for(JZ=JZMIN; JZ<=JZMAX; JZ++)
  {
      for(JY=JYMIN; JY<=JYMAX; JY++)
      {
          for(JX=JXMIN; JX<=JXMAX; JX++)
          {

				/**************** restringiendo rango en (y,z) ****************/ 
                //f1 = Rbis/d;                    // f1(y): y = R;
				f1 = round(Rbis/d);                    // f1(y): y = R;
				f2 = round(sqrt(3)*Rbis/d);            // f2(z): z = sqrt(3)*R
				/**************************************************************/
           
				/***** definiendo el interior los puntos dentro de las esferas *****/

                NP1 = pow(JX-Nanoparticle1.X(), 2) + pow(JY-Nanoparticle1.Y(), 2) + pow(JZ-Nanoparticle1.Z(), 2);   // NP 1
                NP2 = pow(JX-Nanoparticle2.X(), 2) + pow(JY-Nanoparticle2.Y(), 2) + pow(JZ-Nanoparticle2.Z(), 2);   // NP 2
                NP3 = pow(JX-Nanoparticle3.X(), 2) + pow(JY-Nanoparticle3.Y(), 2) + pow(JZ-Nanoparticle3.Z(), 2);   // NP 3
                NP4 = pow(JX-Nanoparticle4.X(), 2) + pow(JY-Nanoparticle4.Y(), 2) + pow(JZ-Nanoparticle4.Z(), 2);   // NP 4
                NP5 = pow(JX-Nanoparticle5.X(), 2) + pow(JY-Nanoparticle5.Y(), 2) + pow(JZ-Nanoparticle5.Z(), 2);   // NP 5
                NP6 = pow(JX-Nanoparticle6.X(), 2) + pow(JY-Nanoparticle6.Y(), 2) + pow(JZ-Nanoparticle6.Z(), 2);   // NP 6

           
				// if(  (fabs(JY) < f1)  &&  (fabs(JZ) < f2)  )                        // version 1 
				 if(   (-f1 < JY) && (JY <= f1)    &&    (-f2 < JZ) && (JZ < f2)   )    // version 2
                // if(  (fabs(JY) < f1)  &&  (fabs(JZ) <= f2)  )                        // version 3 
				{

		                  /***** condiciones de que los puntos esten dentro de las esferas *****/
                          // nanopartículas  
                          if( (rdipSquare - NP1) > 0 )  // diferencia de cuadrados (diferencia de flotantes)
				          {
				                    JX_data[J]=JX;  JY_data[J]=JY;  JZ_data[J]=JZ;
				                    ICOMP[J]=1;
				                    ++J; 
				                    Nanoparticle1.add_dipole();
			              }

                          if( (rdipSquare - NP2) > 0 )  // diferencia de cuadrados (diferencia de flotantes)
				          {
				                    JX_data[J]=JX;  JY_data[J]=JY;  JZ_data[J]=JZ;
				                    ICOMP[J]=1;
				                    ++J; 
				                    Nanoparticle2.add_dipole();
			              }

                          if( (rdipSquare - NP3) > 0 )  // diferencia de cuadrados (diferencia de flotantes)
				          {
				                    JX_data[J]=JX;  JY_data[J]=JY;  JZ_data[J]=JZ;
				                    ICOMP[J]=1;
				                    ++J; 
				                    Nanoparticle3.add_dipole();
			              }

                          if( (rdipSquare - NP4) > 0 )  // diferencia de cuadrados (diferencia de flotantes)
				          {
				                    JX_data[J]=JX;  JY_data[J]=JY;  JZ_data[J]=JZ;
				                    ICOMP[J]=1;
				                    ++J; 
				                    Nanoparticle4.add_dipole();
			              }
                          
                          if( (rdipSquare - NP5) > 0 )  // diferencia de cuadrados (diferencia de flotantes)
				          {
				                    JX_data[J]=JX;  JY_data[J]=JY;  JZ_data[J]=JZ;
				                    ICOMP[J]=1;
				                    ++J; 
				                    Nanoparticle5.add_dipole();
			              }

                          if( (rdipSquare - NP6) > 0 )  // diferencia de cuadrados (diferencia de flotantes)
				          {
				                    JX_data[J]=JX;  JY_data[J]=JY;  JZ_data[J]=JZ;
				                    ICOMP[J]=1;
				                    ++J; 
				                    Nanoparticle6.add_dipole();
			              }

          		}

          }

      }

  }
 /*************************************************************************************/
  /*************************************************************************************/

  /**********************************/
  Vtot = J*pow(d,3);  // N*d^3
  aeff_dip =  cbrt( (3*Vtot)/(4*pi) );  // radio efectivo
  /**********************************/

  /***** Para determinar el rango "real" de los dipolos que se ocuparon en la geometría *****/
	  XDmin = *std::min_element(JX_data.begin(), JX_data.end());  // el * es porque la función regresa la dirección del máximo elemento
	  XDmax = *std::max_element(JX_data.begin(), JX_data.end());

	  YDmin = *std::min_element(JY_data.begin(), JY_data.end());
	  YDmax = *std::max_element(JY_data.begin(), JY_data.end());
	  
	  ZDmin = *std::min_element(JZ_data.begin(), JZ_data.end());
	  ZDmax = *std::max_element(JZ_data.begin(), JZ_data.end());  
  /*****************************************************************************************/

  printf("===================================================================\n");
  printf("\t VALORES CALCULADOS\n");
  printf("-------------------------------------------------------------------\n");
  printf("Calculated dipoles for NP 1:              nNP1 = %d \n", Nanoparticle1.Ndipoles());
  printf("Calculated dipoles for NP 2:              nNP2 = %d \n", Nanoparticle2.Ndipoles());
  printf("Calculated dipoles for NP 3:              nNP3 = %d \n", Nanoparticle3.Ndipoles());
  printf("Calculated dipoles for NP 4:              nNP4 = %d \n", Nanoparticle4.Ndipoles());
  printf("Calculated dipoles for NP 5:              nNP5 = %d \n", Nanoparticle5.Ndipoles());
  printf("Calculated dipoles for NP 6:              nNP6 = %d \n", Nanoparticle6.Ndipoles());
  printf("Volume NP 1                               vNP1 = %1.3E nm^3\n", Nanoparticle1.Ndipoles()*pow(d,3)); 
  printf("Volume NP 2                               vNP2 = %1.3E nm^3\n", Nanoparticle2.Ndipoles()*pow(d,3)); 
  printf("Volume NP 3                               vNP3 = %1.3E nm^3\n", Nanoparticle3.Ndipoles()*pow(d,3)); 
  printf("Volume NP 4                               vNP4 = %1.3E nm^3\n", Nanoparticle4.Ndipoles()*pow(d,3)); 
  printf("Volume NP 5                               vNP5 = %1.3E nm^3\n", Nanoparticle5.Ndipoles()*pow(d,3)); 
  printf("Volume NP 6                               vNP6 = %1.3E nm^3\n\n", Nanoparticle6.Ndipoles()*pow(d,3)); 

  printf("Número total de dipolos calculados:         N  = %d \n", J);
 // printf("Volume of one sphere:                     Vtot = %1.3E nm^3\n", Vtot);
  printf("                                (XDmin, XDmax) = (%d, %d) \n", XDmin, XDmax);
  printf("                                (YDmin, YDmax) = (%d, %d) \n", YDmin, YDmax);
  printf("                                (ZDmin, ZDmax) = (%d, %d) \n", ZDmin, ZDmax);

  printf("-------------------------------------------------------------------\n");
  printf("PARAMETROS PARA DDSCAT\n");
  printf("-------------------------------------------------------------------\n");
  printf("Efective radius:                        aeff = %.2f nm = %.5f um\n", aeff_dip, aeff_dip/1000);
  printf("Periodicidad en y:                        Py = %.2f nm   Py/d = %.2f -> %.f\n", Py, Py/d, round(Py/d));
  printf("Periodicidad en z:                        Pz = %.2f nm   Pz/d = %.2f -> %.0f\n", Pz, Pz/d, round(Pz/d));
  printf("(NO es necesario que Pz = |Zmax-Zmin| porque no se truncan las 2 esferas de arriba y abajo) \n");
  printf("Estimación de la memoria mínima:        MMAX = %.0f      MMAXY = %.0f     MMAXZ = %.0f  \n", 3*Rbis/d, 3*Rbis/d, 4*Rbis/d);
  printf("===================================================================\n");
   printf("Output file: shape.dat\n\n");
   /************************************************************************************/

  /******************************** file output **************************************/
  fprintf(pFile, "---- monolayer of nps: dNP= %.2f nm, Dvirt= %.2f nm, d= %.2f nm ---- \n", dNP, D, d);
  fprintf(pFile, "              %i    = NAT \n", J);
  fprintf(pFile, "1.000  0.000  0.000 = taget vector a1 (in TF)\n");
  fprintf(pFile, "0.000  1.000  0.000 = taget vector a2 (in TF)\n");
  fprintf(pFile, "1.     1.     1.    = lattice spacings (d_x, d_y, d_z)/d (normally 1 1 1) \n" );
  fprintf(pFile, "%.2f   %.2f   %.2f  = x0(1-3) = location in lattice of \"target origin\" \n", X0/d, Y0/d, Z0/d);
   // fprintf(pFile, "0.     0.     0.    = x0(1-3) = location in lattice of \"target origin\"\n");
  fprintf(pFile, "      J    JX      JY      JZ     ICOMP(x,y,z) \n");

  for(int i=0; i<J; i++)
  {
   
      fprintf(pFile, "%7i\t %4i\t %4i\t %4i\t %2i\t %2i\t %2i\n", i+1, JX_data[i], JY_data[i], JZ_data[i], ICOMP[i], ICOMP[i], ICOMP[i]);
  }


   fclose (pFile);

   printf("Fin del programa\n");

  return 0;
}
