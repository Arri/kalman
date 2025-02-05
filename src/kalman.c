/*
Kalman Filter algorithm...
*/
//#include <stdio.h>
#include "main.h"

//--------------------------------------------------------
void initKalman( void )
{
/*
 * Initialization of the kalman filter matrices and vectors:
 */
  printf("[Info] Initializing ...\n");
}

//--------------------------------------------------------

void predict( float u, float v , float* pr_x, float* pr_y,
            float* F, float* B, float* P, float* Q )
{
/*
 * Kalman prediction function
 */
  // Inittializing for now !!!!
  u = 0.0;
  v = 0.0;
  //printf("In Predict \n ");

  *pr_x = *addV2V( dotM2V( F, pr_x ), dotV2S( B, u ) );
  *pr_y = *addV2V( dotM2V( F, pr_y ), dotV2S( B, v ) );
  *P = *addM2M( dotM2M( dotM2M( F, P ), MatT( F ) ), Q );
  return;
}

//
//---------------------------------------------------------

void update( float zx, float zy, float* pr_x, float* pr_y,
            float* P, float* H, float* R, float* I )
{
  /*
   * Updating the kalman parameters based on new measurement values (zx, zy)
   */
  float yx = zx - dotV2V_S( H, pr_x );
  float yy = zy - dotV2V_S( H, pr_y );

  float S  = *R + dotV2V_S( H, dotM2V( P, H ) );  
  float *K = dotV2S( dotM2V( P, H ), (1.0/S) );

  *pr_x = *addV2V( pr_x, dotV2S( K, yx ));
  *pr_y = *addV2V( pr_y, dotV2S( K, yy ));

  *P = *dotM2M( 
             dotM2M( 
                    subsM2M( I, dotV2V_M( K, H ) )
                    , P ), 
             addM2M( 
                    MatT( 
                         subsM2M( I, dotV2V_M( K, H )))
                    , dotV2V_M(
                               dotV2S( K, *R ), K))
                           );

  return;
}
