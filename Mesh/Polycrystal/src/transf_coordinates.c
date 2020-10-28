#include <math.h>
#include <stdio.h>

float matmul(float *omega,float *vec, float *vec2)
{
    int i, j;
    float val;

    for (i=0; i<3; i++){
        val = 0;
        for (j=0; j<3; j++){
            val = val + omega[i*3+j]*vec[j];
        }
        vec2[i] = val;
    }

}

float transform_vector(float *normal,float *vector,float *vector_tt)
{

    float nx, ny, nz, theta_x, theta_y, theta_z, m, n;
    float vector_t[3], normal_t[3], omega_x[9], omega_y[9], omega_z[9];
    float pi=3.141592653589793;
    
    // components of normal vector
    nx = normal[0]; ny = normal[1]; nz = normal[2];

    // Rotation matrix in x
    omega_x[0] = 1; omega_x[1] = 0; omega_x[2] = 0;
    omega_x[3] = 0; omega_x[4] = 0; omega_x[5] = 0;
    omega_x[6] = 0; omega_x[7] = 0; omega_x[8] = 0;

    // Rotation matrix in y
    omega_y[0] = 0; omega_y[1] = 0; omega_y[2] = 0;
    omega_y[3] = 0; omega_y[4] = 1; omega_y[5] = 0;
    omega_y[6] = 0; omega_y[7] = 0; omega_y[8] = 0;

    // Rotation matrix in z
    omega_z[0] = 0; omega_z[1] = 0; omega_z[2] = 0;
    omega_z[3] = 0; omega_z[4] = 0; omega_z[5] = 0;
    omega_z[6] = 0; omega_z[7] = 0; omega_z[8] = 1;

    // (I) and (VII) 
    if (nx>0 && ny>0 && nz>0 || nx<0 && ny<0 && nz<0){
        // Rotation in z
        theta_z = atan(ny/nx);
        m = cos(theta_z); n = sin(theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,normal,normal_t);
        matmul(omega_z,vector,vector_t); 
        // Rotation in y
        theta_y = atan(normal_t[0]/normal_t[2]);
        m = cos(theta_y); n = sin(theta_y);
        omega_y[0] = m; omega_y[2] = -n; 
        omega_y[6] = n; omega_y[8] =  m;
        matmul(omega_y,vector_t,vector_tt); 
    }

    // (II) and (VIII)
    if (nx<0 && ny>0 && nz>0 || nx>0 && ny<0 && nz<0){
        // Rotation in z
        theta_z = atan(-nx/ny);
        m = cos(theta_z); n = sin(theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,normal,normal_t); 
        matmul(omega_z,vector,vector_t); 
        // Rotation in x
        theta_x = atan(normal_t[1]/normal_t[2]);
        m = cos(-theta_x); n = sin(-theta_x);
        omega_x[4] =  m; omega_x[5] = n; 
        omega_x[7] = -n; omega_x[8] = m;
        matmul(omega_x,vector_t,vector_tt);
        
    }

    // (III) and (V)
    if (nx<0 && ny<0 && nz>0 || nx>0 && ny>0 && nz<0){
        // Rotation in z
        theta_z = atan(ny/nx);
        m = cos(theta_z); n = sin(theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,normal,normal_t);
        matmul(omega_z,vector,vector_t); 
        // Rotation in y
        theta_y = atan(-normal_t[0]/normal_t[2]);
        m = cos(-theta_y); n = sin(-theta_y);
        omega_y[0] = m; omega_y[2] = -n; 
        omega_y[6] = n; omega_y[8] =  m;
        matmul(omega_y,vector_t,vector_tt); 
    }

    // (IV) and (VI)
    if (nx>0 && ny<0 && nz>0 || nx<0 && ny>0 && nz<0){
        // Rotation in z
        theta_z = atan(-nx/ny);
        m = cos(theta_z); n = sin(theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,normal,normal_t); 
        matmul(omega_z,vector,vector_t); 
        // Rotation in x
        theta_x = atan(-normal_t[1]/normal_t[2]);
        m = cos(theta_x); n = sin(theta_x);
        omega_x[4] =  m; omega_x[5] = n; 
        omega_x[7] = -n; omega_x[8] = m;
        matmul(omega_x,vector_t,vector_tt);
    }

    // (Y,Z)
    if (nx==0){
        // (I) and (III)
        if (ny>0 && nz>0 || ny<0 && nz<0){
            // Rotation in x
            theta_x = atan(ny/nz);
            m = cos(-theta_x); n = sin(-theta_x);
            omega_x[4] =  m; omega_x[5] = n; 
            omega_x[7] = -n; omega_x[8] = m; 
            matmul(omega_x,vector,vector_tt);  
        }
        // (II) and (IV)
        if (ny<0 && nz>0 || ny>0 && nz<0){
            // Rotation in x
            theta_x = atan(-ny/nz);
            m = cos(theta_x); n = sin(theta_x);
            omega_x[4] =  m; omega_x[5] = n; 
            omega_x[7] = -n; omega_x[8] = m; 
            matmul(omega_x,vector,vector_tt); 
        }
    }

    // (X,Z)
    if (ny==0){
        // (I) and (III)
        if (nx>0 && nz>0 || nx<0 && nz<0){
            // Rotation in y
            theta_y = atan(nx/nz);
            m = cos(theta_y); n = sin(theta_y);
            omega_y[0] = m; omega_y[2] = -n; 
            omega_y[6] = n; omega_y[8] =  m;
            matmul(omega_y,vector,vector_tt);
        }
        // (II) and (IV)
        if (nx<0 && nz>0 || nx>0 && nz<0){
            // Rotation in y
            theta_y = atan(-nx/nz);
            m = cos(-theta_y); n = sin(-theta_y);
            omega_y[0] = m; omega_y[2] = -n; 
            omega_y[6] = n; omega_y[8] =  m;
            matmul(omega_y,vector,vector_tt);  
        }
    }

    // (X,Y)
    if (nz==0){
        // (I) and (III)
        if (nx>0 && ny>0 || nx<0 && ny<0){
            // Rotation in z
            theta_z = atan(ny/nx);
            m = cos(theta_z); n = sin(theta_z);
            omega_z[0] =  m; omega_z[1] = n; 
            omega_z[3] = -n; omega_z[4] = m; 
            matmul(omega_z,vector,vector_t); 
            // Rotation in y
            theta_y = pi/2;
            m = cos(theta_y); n = sin(theta_y);
            omega_y[0] = m; omega_y[2] = -n; 
            omega_y[6] = n; omega_y[8] =  m;
            matmul(omega_y,vector_t,vector_tt);
        }
        // (II) and (IV)
        if (nx<0 && ny>0 || nx>0 && ny<0){
            // Rotation in z
            theta_z = atan(-nx/ny);
            m = cos(theta_z); n = sin(theta_z);
            omega_z[0] =  m; omega_z[1] = n; 
            omega_z[3] = -n; omega_z[4] = m; 
            matmul(omega_z,vector,vector_t); 
            // Rotation in x
            theta_x = pi/2;
            m = cos(-theta_x); n = sin(-theta_x);
            omega_x[4] =  m; omega_x[5] = n; 
            omega_x[7] = -n; omega_x[8] = m;
            matmul(omega_x,vector_t,vector_tt);
        }
    }

    // Axis x
    if (ny==0 && nz==0){
        // Rotation in y
        theta_y = pi/2;
        m = cos(-theta_y); n = sin(-theta_y);
        omega_y[0] = m; omega_y[2] = -n; 
        omega_y[6] = n; omega_y[8] =  m;
        matmul(omega_y,vector,vector_tt);    
    }

    // Axis y
    if (nx==0 && nz==0){
        // Rotation in x
        theta_x = pi/2;
        m = cos(theta_x); n = sin(theta_x);
        omega_x[4] =  m; omega_x[5] = n; 
        omega_x[7] = -n; omega_x[8] = m; 
        matmul(omega_x,vector,vector_tt);

    }

    // Axis z
    if (nx==0 && ny==0){
        vector_tt[0] = vector[0];
        vector_tt[1] = vector[1];
        vector_tt[2] = vector[2];
    }

}
    
float transform_vector_inverse(float *normal,float *vector_tt,float *vector)
{

    float nx, ny, nz, theta_x, theta_y, theta_z, m, n;
    float vector_t[3], normal_t[3], omega_x[9], omega_y[9], omega_z[9];
    float pi=3.141592653589793;
    
    // components of normal vector
    nx = normal[0]; ny = normal[1]; nz = normal[2];

    // Rotation matrix in x
    omega_x[0] = 1; omega_x[1] = 0; omega_x[2] = 0;
    omega_x[3] = 0; omega_x[4] = 0; omega_x[5] = 0;
    omega_x[6] = 0; omega_x[7] = 0; omega_x[8] = 0;

    // Rotation matrix in y
    omega_y[0] = 0; omega_y[1] = 0; omega_y[2] = 0;
    omega_y[3] = 0; omega_y[4] = 1; omega_y[5] = 0;
    omega_y[6] = 0; omega_y[7] = 0; omega_y[8] = 0;

    // Rotation matrix in z
    omega_z[0] = 0; omega_z[1] = 0; omega_z[2] = 0;
    omega_z[3] = 0; omega_z[4] = 0; omega_z[5] = 0;
    omega_z[6] = 0; omega_z[7] = 0; omega_z[8] = 1;

    // (I) and (VII) 
    if (nx>0 && ny>0 && nz>0 || nx<0 && ny<0 && nz<0){
        // Angle in z
        theta_z = atan(ny/nx);
        m = cos(theta_z); n = sin(theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,normal,normal_t); 
        // Rotation in y
        theta_y = atan(normal_t[0]/normal_t[2]);
        m = cos(-theta_y); n = sin(-theta_y);
        omega_y[0] = m; omega_y[2] = -n; 
        omega_y[6] = n; omega_y[8] =  m;
        matmul(omega_y,vector_tt,vector_t); 
        // Rotation in z
        m = cos(-theta_z); n = sin(-theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,vector_t,vector);
    }

    // (II) and (VIII)
    if (nx<0 && ny>0 && nz>0 || nx>0 && ny<0 && nz<0){
        // Angle in z
        theta_z = atan(-nx/ny);
        m = cos(theta_z); n = sin(theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,normal,normal_t);  
        // Rotation in x
        theta_x = atan(normal_t[1]/normal_t[2]);
        m = cos(theta_x); n = sin(theta_x);
        omega_x[4] =  m; omega_x[5] = n; 
        omega_x[7] = -n; omega_x[8] = m;
        matmul(omega_x,vector_tt,vector_t);
        // Rotation in z
        m = cos(-theta_z); n = sin(-theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,vector_t,vector);      
    }

    // (III) and (V)
    if (nx<0 && ny<0 && nz>0 || nx>0 && ny>0 && nz<0){
        // Angle in z
        theta_z = atan(ny/nx);
        m = cos(theta_z); n = sin(theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,normal,normal_t);
        // Rotation in y
        theta_y = atan(-normal_t[0]/normal_t[2]);
        m = cos(theta_y); n = sin(theta_y);
        omega_y[0] = m; omega_y[2] = -n; 
        omega_y[6] = n; omega_y[8] =  m;
        matmul(omega_y,vector_tt,vector_t); 
        // Rotation in z
        m = cos(-theta_z); n = sin(-theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,vector_t,vector);
    }

    // (IV) and (VI)
    if (nx>0 && ny<0 && nz>0 || nx<0 && ny>0 && nz<0){
        // Rotation in z
        theta_z = atan(-nx/ny);
        m = cos(theta_z); n = sin(theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,normal,normal_t); 
        // Rotation in x
        theta_x = atan(-normal_t[1]/normal_t[2]);
        m = cos(-theta_x); n = sin(-theta_x);
        omega_x[4] =  m; omega_x[5] = n; 
        omega_x[7] = -n; omega_x[8] = m;
        matmul(omega_x,vector_tt,vector_t);
        // Rotation in z
        theta_z = atan(-nx/ny);
        m = cos(-theta_z); n = sin(-theta_z);
        omega_z[0] =  m; omega_z[1] = n; 
        omega_z[3] = -n; omega_z[4] = m; 
        matmul(omega_z,vector_t,vector); 
    }

    // (Y,Z)
    if (nx==0){
        // (I) and (III)
        if (ny>0 && nz>0 || ny<0 && nz<0){
            // Rotation in x
            theta_x = atan(ny/nz);
            m = cos(theta_x); n = sin(theta_x);
            omega_x[4] =  m; omega_x[5] = n; 
            omega_x[7] = -n; omega_x[8] = m; 
            matmul(omega_x,vector_tt,vector);  
        }
        // (II) and (IV)
        if (ny<0 && nz>0 || ny>0 && nz<0){
            // Rotation in x
            theta_x = atan(-ny/nz);
            m = cos(-theta_x); n = sin(-theta_x);
            omega_x[4] =  m; omega_x[5] = n; 
            omega_x[7] = -n; omega_x[8] = m; 
            matmul(omega_x,vector_tt,vector); 
        }
    }

    // (X,Z)
    if (ny==0){
        // (I) and (III)
        if (nx>0 && nz>0 || nx<0 && nz<0){
            // Rotation in y
            theta_y = atan(nx/nz);
            m = cos(-theta_y); n = sin(-theta_y);
            omega_y[0] = m; omega_y[2] = -n; 
            omega_y[6] = n; omega_y[8] =  m;
            matmul(omega_y,vector_tt,vector);
        }
        // (II) and (IV)
        if (nx<0 && nz>0 || nx>0 && nz<0){
            // Rotation in y
            theta_y = atan(-nx/nz);
            m = cos(theta_y); n = sin(theta_y);
            omega_y[0] = m; omega_y[2] = -n; 
            omega_y[6] = n; omega_y[8] =  m;
            matmul(omega_y,vector_tt,vector);  
        }
    }

    // (X,Y)
    if (nz==0){
        // (I) and (III)
        if (nx>0 && ny>0 || nx<0 && ny<0){
            // Rotation in y
            theta_y = pi/2;
            m = cos(-theta_y); n = sin(-theta_y);
            omega_y[0] = m; omega_y[2] = -n; 
            omega_y[6] = n; omega_y[8] =  m;
            matmul(omega_y,vector_tt,vector_t);
            // Rotation in z
            theta_z = atan(ny/nx);
            m = cos(-theta_z); n = sin(-theta_z);
            omega_z[0] =  m; omega_z[1] = n; 
            omega_z[3] = -n; omega_z[4] = m; 
            matmul(omega_z,vector_t,vector);     
        }
        // (II) and (IV)
        if (nx<0 && ny>0 || nx>0 && ny<0){
            // Rotation in x
            theta_x = pi/2;
            m = cos(theta_x); n = sin(theta_x);
            omega_x[4] =  m; omega_x[5] = n; 
            omega_x[7] = -n; omega_x[8] = m;
            matmul(omega_x,vector_tt,vector_t);
            // Rotation in z
            theta_z = atan(-nx/ny);
            m = cos(-theta_z); n = sin(-theta_z);
            omega_z[0] =  m; omega_z[1] = n; 
            omega_z[3] = -n; omega_z[4] = m; 
            matmul(omega_z,vector_t,vector); 
        }
    }

    // Axis x
    if (ny==0 && nz==0){
        // Rotation in y
        theta_y = pi/2;
        m = cos(theta_y); n = sin(theta_y);
        omega_y[0] = m; omega_y[2] = -n; 
        omega_y[6] = n; omega_y[8] =  m;
        matmul(omega_y,vector_tt,vector);    
    }

    // Axis y
    if (nx==0 && nz==0){
        // Rotation in x
        theta_x = pi/2;
        m = cos(-theta_x); n = sin(-theta_x);
        omega_x[4] =  m; omega_x[5] = n; 
        omega_x[7] = -n; omega_x[8] = m; 
        matmul(omega_x,vector_tt,vector);

    }

    // Axis z
    if (nx==0 && ny==0){
        vector[0] = vector_tt[0];
        vector[1] = vector_tt[1];
        vector[2] = vector_tt[2];
    }


}
